/*
 * boxcar.cpp
 *
 *  Created on: Apr 30, 2020
 *      Author: ypmen
 */

#include <stdlib.h>
#include <string.h>
#include <cmath>
#ifdef __AVX2__
	#include <immintrin.h>
#endif

#include "boxcar.h"

Boxcar::Boxcar()
{
	counter = 0;
	tsamp = 0.;
	nsamples = 0;
	dms = 0.;
	ddm = 0.;
	ndm = 0;
	maxwn = 0;
	fmin = 0.;
	fmax = 0.;
	mxS = NULL;
	mxwn = NULL;
}

Boxcar::Boxcar(const Boxcar &boxcar)
{
	counter = boxcar.counter;
	tsamp = boxcar.tsamp;
	nsamples = boxcar.nsamples;
	dms = boxcar.dms;
	ddm = boxcar.ddm;
	ndm = boxcar.ndm;
	maxwn = boxcar.maxwn;
	fmin = boxcar.fmin;
	fmax = boxcar.fmax;
	if (boxcar.mxS != NULL)
	{
		mxS = new float [ndm*nsamples];
		memcpy(mxS, boxcar.mxS, sizeof(float)*ndm*nsamples);
	}
	else
	{
		mxS = NULL;
	}

	if (boxcar.mxwn != NULL)
	{
		mxwn = new int [ndm*nsamples];
		memcpy(mxwn, boxcar.mxwn, sizeof(int)*ndm*nsamples);
	}
	else
	{
		mxwn = NULL;
	}
	finterval = boxcar.finterval;
}

Boxcar & Boxcar::operator=(const Boxcar &boxcar)
{
	counter = boxcar.counter;
	tsamp = boxcar.tsamp;
	nsamples = boxcar.nsamples;
	dms = boxcar.dms;
	ddm = boxcar.ddm;
	ndm = boxcar.ndm;
	maxwn = boxcar.maxwn;
	fmin = boxcar.fmin;
	fmax = boxcar.fmax;
	if (boxcar.mxS != NULL)
	{
		mxS = new float [ndm*nsamples];
		memcpy(mxS, boxcar.mxS, sizeof(float)*ndm*nsamples);
	}
	else
	{
		mxS = NULL;
	}

	if (boxcar.mxwn != NULL)
	{
		mxwn = new int [ndm*nsamples];
		memcpy(mxwn, boxcar.mxwn, sizeof(int)*ndm*nsamples);
	}
	else
	{
		mxwn = NULL;
	}
	finterval = boxcar.finterval;

	return *this;
}

Boxcar::~Boxcar()
{
	if (mxS != NULL)
	{
		delete [] mxS;
		mxS = NULL;
	}

	if (mxwn != NULL)
	{
		delete [] mxwn;
		mxwn = NULL;
	}
}

void Boxcar::resize(long int ns, long int nd)
{
	if (ns != nsamples or nd != ndm)
	{
		if (mxS != NULL) delete [] mxS;
		if (mxwn != NULL) delete [] mxwn;

		nsamples = ns;
		ndm = nd;

		mxS = new float [ndm*nsamples];
		mxwn = new int [ndm*nsamples];

		finterval.resize(ndm);
		for (long int j=0; j<ndm; j++)
		{
			finterval[j].resize(nsamples);
		}
	}

	tsamp = 0.;
	dms = 0.;
	ddm = 0.;
	memset(mxS, 0, sizeof(float)*ndm*nsamples);
	memset(mxwn, 0, sizeof(int)*ndm*nsamples);
}

void Boxcar::prepare(RealTime::SubbandDedispersion &dedisp)
{
	resize(dedisp.ndump+dedisp.noverlap, dedisp.ndm);
	dms = dedisp.dms;
	ddm = dedisp.ddm;
	tsamp = dedisp.tsamp;

	fmin = 1e6;
	fmax = 0.;
	for (long int j=0; j<dedisp.nchans; j++)
	{
		fmax = dedisp.frequencies[j]>fmax? dedisp.frequencies[j]:fmax;
		fmin = dedisp.frequencies[j]<fmin? dedisp.frequencies[j]:fmin;
	}
}

bool Boxcar::run(RealTime::SubbandDedispersion &dedisp, vector<int> &vwn, bool iqr)
{
	counter = dedisp.counter;
	if (counter < dedisp.offset-dedisp.noverlap+dedisp.ndump) return false;

	int n = vwn.size();
	maxwn = vwn[n-1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int k=0; k<dedisp.ndm; k++)
	{
		if (!repeater)
			match(k, vwn, dedisp, iqr);
		else
			match2D(k, vwn, dedisp);
	}

	return true;
}

void Boxcar::match(int idm, vector<int> &vwn, RealTime::SubbandDedispersion &dedisp, bool iqr)
{
	vector<float> tim;
	dedisp.get_timdata(tim, idm, true);

	int n = vwn.size();

#ifdef __AVX2__
	float *vS = (float *)_mm_malloc(sizeof(float)*nsamples, 32);
#else
	float *vS = (float *)malloc(sizeof(float)*nsamples);
#endif
	memset(vS, 0, sizeof(float)*nsamples);

#ifdef __AVX2__
	int *vwn_maxS = (int *)_mm_malloc(sizeof(float)*nsamples, 32);
#else
	int *vwn_maxS = (int *)malloc(sizeof(float)*nsamples);
#endif
	memset(vwn_maxS, 0, sizeof(int)*nsamples);

	float mean = 0.;
	float var = 0.;
	
	if (!iqr)
	{
		for (long int i=0; i<nsamples; i++)
		{
			mean += tim[i];
			var += tim[i]*tim[i];
		}
		mean /= nsamples;
		var /= nsamples;
		var -= mean*mean;
	}
	else
	{
		vector<float> timcopy = tim;
		std::nth_element(timcopy.begin(), timcopy.begin()+nsamples/4, timcopy.end(), std::less<float>());
		float Q1 = timcopy[nsamples/4];
		std::nth_element(timcopy.begin(), timcopy.begin()+nsamples/2, timcopy.end(), std::less<float>());
		float Q2 = timcopy[nsamples/2];
		std::nth_element(timcopy.begin(), timcopy.begin()+nsamples/4, timcopy.end(), std::greater<float>());
		float Q3 = timcopy[nsamples/4];

		mean = Q2;
		var = ((Q3-Q1)/1.349)*((Q3-Q1)/1.349);
	}

#ifdef __AVX2__
	float *csump = (float *)_mm_malloc(sizeof(float)*nsamples, 32);
#else
	float *csump = (float *)malloc(sizeof(float)*nsamples);
#endif

	float csum = 0.;
	for (long int i=0; i<nsamples; i++)
	{
		csum += tim[i] - mean;
		csump[i] = csum;
	}

	if (var > 0)
	{
		for (long int k=0; k<n; k++)
		{
			int wn = vwn[k];
			int wl = wn/2+1;
			int wh = (wn-1)/2;
			
			float temp = sqrt(1./(wn*var));	

#ifndef __AVX2__
			for (long int i=wl; i<nsamples-wh; i++)
			{
				float boxsum = csump[i+wh]-csump[i-wl];
				float S = boxsum*temp;	
				if (S>=vS[i])
				{
					vS[i] = S;
					vwn_maxS[i] = wn;
				}
			}
#else
			int len = nsamples-wh-wl;
			__m256 avx_temp = _mm256_set_ps(temp, temp, temp, temp, temp, temp, temp, temp);
			__m256 avx_wn = _mm256_set_ps(wn, wn, wn, wn, wn, wn, wn, wn);
			for (long int i=0; i<len/8; i++)
			{
				__m256 avx_boxsum = _mm256_sub_ps(_mm256_loadu_ps(csump+wl+wh+i*8), _mm256_load_ps(csump+i*8));
				__m256 avx_S = _mm256_mul_ps(avx_boxsum, avx_temp);
				__m256 avx_vS = _mm256_loadu_ps(vS+wl+i*8);
				__m256 avx_vwn_maxS = _mm256_cvtepi32_ps(_mm256_loadu_si256((__m256i *)(vwn_maxS+wl+i*8)));
				__m256 avx_mask = _mm256_cmp_ps(avx_S, avx_vS, 13);
				__m256 avx_maskinv = _mm256_cmp_ps(avx_S, avx_vS, 1);
				avx_vS = _mm256_add_ps(_mm256_and_ps(avx_mask, avx_S), _mm256_and_ps(avx_maskinv, avx_vS));
				_mm256_storeu_ps(vS+wl+i*8, avx_vS);
				avx_vwn_maxS = _mm256_add_ps(_mm256_and_ps(avx_mask, avx_wn), _mm256_and_ps(avx_maskinv, avx_vwn_maxS));
				_mm256_storeu_si256((__m256i *)(vwn_maxS+wl+i*8), _mm256_cvtps_epi32(avx_vwn_maxS));
			}
#endif
		}
	}

	memcpy(mxS+idm*nsamples, vS, sizeof(float)*nsamples);
	memcpy(mxwn+idm*nsamples, vwn_maxS, sizeof(int)*nsamples);

#ifdef __AVX2__
	_mm_free(csump);
	_mm_free(vS);
	_mm_free(vwn_maxS);
#else
	free(csump);
	free(vS);
	free(vwn_maxS);
#endif
}

//not work
void Boxcar::match2D(int idm, vector<int> &vwn, RealTime::SubbandDedispersion &dedisp)
{
	int n = vwn.size();

	float *vS = new float [nsamples];
	memset(vS, 0, sizeof(float)*nsamples);

	int *vwn_maxS = new int [nsamples];
	memset(vwn_maxS, 0, sizeof(int)*nsamples);

	int *vfl = new int [nsamples];
	int *vfh = new int [nsamples];

	double *submean = new double [dedisp.sub.nchans];
	double *subvar = new double [dedisp.sub.nchans];
	memset(submean, 0, sizeof(double)*dedisp.sub.nchans);
	memset(subvar, 0, sizeof(double)*dedisp.sub.nchans);

	vector<float> sub;
	dedisp.get_subdata(sub, idm);
	float *psd = &sub[0];
	for (long int i=0; i<nsamples; i++)
	{
		vfl[i] = 0;
		vfh[i] = dedisp.sub.nchans;
		for (long int j=0; j<dedisp.sub.nchans; j++)
		{
			submean[j] += *psd;
			subvar[j] += (*psd)*(*psd);
			psd++;
		}
	}

	for (long int j=0; j<dedisp.sub.nchans; j++)
	{
		submean[j] /= nsamples;
		subvar[j] /= nsamples;
		subvar[j] -= submean[j]*submean[j];
	}

	float *mxp = new float [nsamples*dedisp.sub.nchans];
	for (long int i=0; i<nsamples; i++)
	{
		for (long int j=0; j<dedisp.sub.nchans; j++)
		{
			*mxp -= submean[j];
			mxp++;
		}
	}

	float *fcum = new float [dedisp.sub.nchans];

	for (long int k=0; k<n; k++)
	{
		int wn = vwn[k];
		int wl = wn/2+1;
		int wh = (wn-1)/2;
		int fl = -1;
		int fh = dedisp.sub.nchans;
		double boxsum = 0.;

		memset(fcum, 0, sizeof(float)*dedisp.sub.nchans);
		psd = mxp;
		long int m=0;
		for (long int l=0; l<wh; l++)
		{
			float cum=0.;
			for (long int j=0; j<dedisp.sub.nchans; j++)
			{
				cum += psd[m++];
				fcum[j] = cum;
			}
		}

		for (long int i=0; i<nsamples; i++)
		{
			if (i-wl >= 0 and i+wh < nsamples)
			{
				float cuml=0.;
				float cumh=0.;
				float *pl = psd-wl*dedisp.sub.nchans;
				float *ph = psd+wh*dedisp.sub.nchans;
				for (long int j=0; j<dedisp.sub.nchans; j++)
				{
					cuml += pl[j];
					cumh += ph[j];
					fcum[j] -= cuml;
					fcum[j] += cumh;
				}
			}
			else if (i-wl < 0 and i+wh < nsamples)
			{
				float cumh=0.;
				float *ph = psd+wh*dedisp.sub.nchans;
				for (long int j=0; j<dedisp.sub.nchans; j++)
				{
					cumh += ph[j];
					fcum[j] += cumh;
				}
			}
			else if (i+wh >= nsamples and i-wl >= 0)
			{
				float cuml=0.;
				float *pl = psd-wl*dedisp.sub.nchans;
				for (long int j=0; j<dedisp.sub.nchans; j++)
				{
					cuml += pl[j];
					fcum[j] -= cuml;
				}				
			}
			else
			{
				cerr<<"Error: width is larger than length"<<endl;
				exit(-1);
			}

			float var = 0;
			float maxfcum = fcum[dedisp.sub.nchans-1];
			for (long int j=0; j<dedisp.sub.nchans; j++)
			{
				float level = fcum[j]/maxfcum;

				if (level >= LOWLEVEL and fl==-1)
				{
					fl = j;
				}

				if (level > HIGHLEVEL and fh==dedisp.sub.nchans)
				{
					fh = j;
				}

				if (level >= LOWLEVEL and level <= HIGHLEVEL)
				{
					var += subvar[j];
				}
			}

			if (fl>0)
				boxsum = fcum[fh-1]-fcum[fl-1];
			else
				boxsum = fcum[fh-1];

			float S = boxsum*boxsum/(wn*var);
			if (S>=vS[i])
			{
				vS[i] = S;
				vwn_maxS[i] = wn;
				vfl[i] = fl;
				vfh[i] = fh;
			}
			psd += dedisp.sub.nchans;
		}
	}

	delete [] fcum;

	memcpy(mxS+idm*nsamples, vS, sizeof(float)*nsamples);
	memcpy(mxwn+idm*nsamples, vwn_maxS, sizeof(int)*nsamples);

	for (long int i=0; i++; i<nsamples)
	{
		finterval[idm][i] = pair<int, int>(vfl[i], vfh[i]);
	}

	delete [] mxp;
	delete [] submean;
	delete [] subvar;
	delete [] vfh;
	delete [] vfl;
	delete [] vwn_maxS;
	delete [] vS;
}

void Boxcar::dump2txt(const string fname) const
{
	FILE *fptr = fopen(fname.c_str(), "w");

	float *pS = mxS;
	for (long int j=0; j<ndm; j++)
	{
		for (long int i=0; i<nsamples; i++)
		{
			fprintf(fptr,"%f ", *pS);
			pS++;
		}
		fprintf(fptr,"\n");
	}

	fclose(fptr);
}
