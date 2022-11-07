/*
 * rfi.cpp
 *
 *  Created on: May 6, 2020
 *      Author: ypmen
 */

#include "string.h"

#ifdef __AVX2__
#include "avx2.h"
#endif

#include "rfi.h"
#include "kdtree.h"
#include "dedisperse.h"
#include "logging.h"
#include <random>

using namespace std;

RFI::RFI()
{
	filltype = "mean";
}

RFI::RFI(const RFI &rfi) : DataBuffer<float>(rfi)
{
	weights = rfi.weights;
	filltype = rfi.filltype;
}

RFI & RFI::operator=(const RFI &rfi)
{
	DataBuffer<float>::operator=(rfi);

	weights = rfi.weights;
	filltype = rfi.filltype;

	return *this;  
}

RFI::~RFI(){}

void RFI::prepare(DataBuffer<float> &databuffer)
{
	equalize.prepare(databuffer);
	equalize.close();

	nsamples = databuffer.nsamples;
	nchans = databuffer.nchans;

	resize(nsamples, nchans);

	tsamp = databuffer.tsamp;
	frequencies = databuffer.frequencies;

	weights.resize(nchans, 1);
}

DataBuffer<float> * RFI::zap(DataBuffer<float> &databuffer, const vector<pair<double, double>> &zaplist)
{
	if (zaplist.empty())
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"zapping channels";

	fill(weights.begin(), weights.end(), 1);

	for (long int j=0; j<nchans; j++)
	{
		for (auto k=zaplist.begin(); k!=zaplist.end(); ++k)
		{
			if (frequencies[j]>=(*k).first and frequencies[j]<=(*k).second)
			{
				weights[j] = 0;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples; i++)
	{   
		for (long int j=0; j<nchans; j++)
		{
			databuffer.buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]*weights[j];
		}
	}

	counter += nsamples;

	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return databuffer.get();
}

DataBuffer<float> * RFI::zdot(DataBuffer<float> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform zero-dm matched filter";

#ifndef __AVX2__
	vector<double> xe(nchans, 0.);
	vector<double> xs(nchans, 0.);
	vector<double> alpha(nchans, 0.);
	vector<double> beta(nchans, 0.);
	vector<double> s(nsamples, 0.);
	double se = 0.;
	double ss = 0.;

	for (long int i=0; i<nsamples; i++)
	{
		double temp = 0.;
		for (long int j=0; j<nchans; j++)
		{
			temp += databuffer.buffer[i*nchans+j];
		}

		for (long int j=0; j<nchans; j++)
		{
			xe[j] += databuffer.buffer[i*nchans+j];
		}

		temp /= nchans;
		se += temp;
		ss += temp*temp;
		
		for (long int j=0; j<nchans; j++)
		{
			xs[j] += databuffer.buffer[i*nchans+j]*temp;
		}

		s[i] = temp;
	}
#else
	vector<double, boost::alignment::aligned_allocator<double, 32>> xe(nchans, 0.);
	vector<double, boost::alignment::aligned_allocator<double, 32>> xs(nchans, 0.);
	vector<float, boost::alignment::aligned_allocator<float, 32>> alpha(nchans, 0.);
	vector<float, boost::alignment::aligned_allocator<float, 32>> beta(nchans, 0.);
	vector<float, boost::alignment::aligned_allocator<float, 32>> s(nsamples, 0.);
	double se = 0.;
	double ss = 0.;

	if (nchans % 4 == 0)
	{
		for (long int i=0; i<nsamples; i++)
		{
			double temp = PulsarX::reduce(databuffer.buffer.data()+i*nchans, nchans);

			temp /= nchans;

			PulsarX::accumulate_mean(xe.data(), xs.data(), temp, databuffer.buffer.data()+i*nchans, nchans);

			se += temp;
			ss += temp*temp;

			s[i] = temp;
		}
	}
	else
	{
		for (long int i=0; i<nsamples; i++)
		{
			double temp = 0.;
			for (long int j=0; j<nchans; j++)
			{
				temp += databuffer.buffer[i*nchans+j];
			}

			for (long int j=0; j<nchans; j++)
			{
				xe[j] += databuffer.buffer[i*nchans+j];
			}

			temp /= nchans;
			se += temp;
			ss += temp*temp;
			
			for (long int j=0; j<nchans; j++)
			{
				xs[j] += databuffer.buffer[i*nchans+j]*temp;
			}

			s[i] = temp;
		}
	}
#endif

	double tmp = se*se-ss*nsamples;
	if (tmp != 0)
	{
		for (long int j=0; j<nchans; j++)
		{
			alpha[j] = (xe[j]*se-xs[j]*nsamples)/tmp;
			beta[j] = (xs[j]*se-xe[j]*ss)/tmp;
		}
	}

#ifndef __AVX2__

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			databuffer.buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];
		}
	}

#else
	if (nchans % 8 == 0)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::remove_baseline(databuffer.buffer.data()+i*nchans, databuffer.buffer.data()+i*nchans, alpha.data(), beta.data(), s[i], nchans);
		}
	}
	else
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];
			}
		}
	}
#endif

	databuffer.equalized = false;
	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return databuffer.get();
}

DataBuffer<float> * RFI::zero(DataBuffer<float> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform zero-dm filter";

	if (closable) open();

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples; i++)
	{
		double s = 0;
		for (long int j=0; j<nchans; j++)
		{
			s += databuffer.buffer[i*nchans+j];
		}
		s /= nchans;

		for (long int j=0; j<nchans; j++)
		{
			buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-s;
		}
	}

	equalized = false;

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return this;
}

DataBuffer<float> * RFI::mask(DataBuffer<float> &databuffer, float threRFI2, int td, int fd)
{
	if (closable) open();

	long int nsamples_ds = nsamples/td;
	long int nchans_ds = nchans/fd;

	vector<float> buffer_ds(nsamples_ds*nchans_ds, 0.);

	//downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples_ds; i++)
	{
		for (long int n=0; n<td; n++)
		{
			for (long int k=0; k<fd; k++)
			{
				for (long int j=0; j<nchans_ds; j++)
				{
					buffer_ds[i*nchans_ds+j] += databuffer.buffer[(i*td+n)*nchans+j*fd+k];
				}
			}
		}
	}

	buffer = databuffer.buffer;

	vector<float> buffer_dscopy = buffer_ds;
	std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/4, buffer_dscopy.end(), std::less<float>());
	float Q1 = buffer_dscopy[nsamples_ds*nchans_ds/4];
	std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/2, buffer_dscopy.end(), std::less<float>());
	float Q2 = buffer_dscopy[nsamples_ds*nchans_ds/2];
	std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/4, buffer_dscopy.end(), std::greater<float>());
	float Q3 = buffer_dscopy[nsamples_ds*nchans_ds/4];

	float mean = Q2;
	float var = ((Q3-Q1)/1.349)*((Q3-Q1)/1.349);
	float thre = threRFI2*var;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples_ds; i++)
	{
		for (long int n=0; n<td; n++)
		{
			for (long int j=0; j<nchans_ds; j++)
			{
				for (long int k=0; k<fd; k++)
				{
					if ((buffer_ds[i*nchans_ds+j]-mean)*(buffer_ds[i*nchans_ds+j]-mean)>thre)
						buffer[(i*td+n)*nchans+j*fd+k] = mean;
				}
			}
		}
	}

	equalized = databuffer.equalized;

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	return this;
}

DataBuffer<float> * RFI::kadaneF(DataBuffer<float> &databuffer, float threRFI2, double widthlimit, int td, int fd)
{
	equalize.filter(databuffer);

#ifndef __AVX2__
	if (!databuffer.equalized)
	{
		BOOST_LOG_TRIVIAL(warning)<<"Warning: data is not equalize, kadaneF filter will not be performed"<<endl;
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform kadane filter";

	if (closable) open();

	vector<float> bufferT(nchans*nsamples, 0.);

	transpose_pad<float>(&bufferT[0], &databuffer.buffer[0], nsamples, nchans);

	long int nsamples_ds = nsamples/td;
	long int nchans_ds = nchans/fd;

	vector<float> bufferT_ds(nchans_ds*nsamples_ds, 0.);

	//downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchans_ds; j++)
	{
		for (long int k=0; k<fd; k++)
		{
			for (long int n=0; n<td; n++)
			{
				for (long int i=0; i<nsamples_ds; i++)       
				{
					bufferT_ds[j*nsamples_ds+i] += bufferT[(j*fd+k)*nsamples+(i*td+n)];
				}
			}
		}
	}

#ifdef _OPENMP
	float *chdata_t = new float [num_threads*nsamples_ds];
	memset(chdata_t, 0, sizeof(float)*num_threads*nsamples_ds);

	std::vector<std::random_device> r(num_threads);
	std::vector<std::mt19937> generators;
	for (long int k=0; k<num_threads; k++) generators.emplace_back(std::mt19937(r[k]()));
	std::vector<std::normal_distribution<float>> distributions(num_threads, std::normal_distribution<float>(0., 1.));
#else
	float *chdata_t = new float [nsamples_ds];
	memset(chdata_t, 0, sizeof(float)*nsamples_ds);

	std::random_device r;
	std::mt19937 generator(r());
	std::normal_distribution<float> distribution(0., 1.);
#endif

	int wnlimit = widthlimit/tsamp/td;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchans_ds; j++)
	{
#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
		float *chdata = chdata_t+thread_id*nsamples_ds;
#else
		float *chdata = chdata_t;
#endif

		for (long int i=0; i<nsamples_ds; i++)
		{
			chdata[i] = bufferT_ds[j*nsamples_ds+i];
		}

		long int start, end;
		float snr2=0.;
		float boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
		long int wn = end-start+1;

		float mean = 0.;
		float var = td*fd;    
		if (wn > wnlimit)
		{
			snr2 = boxsum*boxsum/(wn*var);
		}
		else
		{
			snr2 = 0.;
		}

		start *= td;
		end += 1;
		end *= td;
		if (snr2 > threRFI2)
		{
			if (filltype == "mean")
			{
				for (long int k=0; k<fd; k++)
				{
					for (long int i=start; i<end; i++)
					{
						bufferT[(j*fd+k)*nsamples+i] = 0.;
					}
				}
			}
			else
			{
				for (long int k=0; k<fd; k++)
				{
					for (long int i=start; i<end; i++)
					{
#ifdef _OPENMP
						bufferT[(j*fd+k)*nsamples+i] = distributions[thread_id](generators[thread_id]);
#else
						bufferT[(j*fd+k)*nsamples+i] = distribution(generator);
#endif
					}
				}
			}
		}

		//<0
		for (long int i=0; i<nsamples_ds; i++)
		{
			chdata[i] = -chdata[i];
		}

		snr2=0.;
		boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
		wn = end-start+1;

		if (wn > wnlimit)
		{
			snr2 = boxsum*boxsum/(wn*var);
		}
		else
		{
			snr2 = 0.;
		}

		start *= td;
		end += 1;
		end *= td;
		if (snr2 > threRFI2)
		{
			if (filltype == "mean")
			{
				for (long int k=0; k<fd; k++)
				{
					for (long int i=start; i<end; i++)
					{
						bufferT[(j*fd+k)*nsamples+i] = 0.;
					}
				}
			}
			else
			{
				for (long int k=0; k<fd; k++)
				{
					for (long int i=start; i<end; i++)
					{
#ifdef _OPENMP
						bufferT[(j*fd+k)*nsamples+i] = distributions[thread_id](generators[thread_id]);
#else
						bufferT[(j*fd+k)*nsamples+i] = distribution(generator);
#endif
					}
				}
			}
		}
		
	}

	transpose_pad<float>(&buffer[0], &bufferT[0], nchans, nsamples);

	equalized = databuffer.equalized;

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	delete [] chdata_t;

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return this;
#else
	if (nchans % 8 == 0)
	{
		if (!databuffer.equalized)
		{
			BOOST_LOG_TRIVIAL(warning)<<"Warning: data is not equalize, kadaneF filter will not be performed"<<endl;
			return databuffer.get();
		}

		BOOST_LOG_TRIVIAL(debug)<<"perform kadane filter";

		int wnlimit = widthlimit/tsamp/td;

		long int nsamples_ds = nsamples/td;
		long int nchans_ds = nchans/fd;

		std::vector<float, boost::alignment::aligned_allocator<float, 32>> buffer_ds(nchans_ds*nsamples_ds, 0.);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples_ds; i++)
		{
			for (long int n=0; n<td; n++)
			{
				for (long int k=0; k<fd; k++)
				{
					for (long int j=0; j<nchans_ds; j++)
					{
						buffer_ds[i*nchans_ds+j] += databuffer.buffer[(i*td+n)*databuffer.nchans+j*fd+k];
					}
				}
			}
		}

		std::vector<float, boost::alignment::aligned_allocator<float, 32>> boxsum(nchans_ds, 0.), snr2(nchans_ds, 0.);
		std::vector<int, boost::alignment::aligned_allocator<int, 32>> start(nchans_ds, 0), end(nchans_ds, 0);
		std::vector<int, boost::alignment::aligned_allocator<int, 32>> wn(nchans_ds, 0);

		PulsarX::kadane2D(boxsum.data(), start.data(), end.data(), buffer_ds.data(), nsamples_ds, nchans_ds);

		for (long int j=0; j<nchans_ds; j++)
		{
			wn[j] = end[j]-start[j]+1;

			float mean = 0.;
			float var = td*fd;

			if (wn[j] > wnlimit)
			{
				snr2[j] = boxsum[j]*boxsum[j]/(wn[j]*var);
			}
			else
			{
				snr2[j] = 0.;
			}

			start[j] *= td;
			end[j] += 1;
			end[j] *= td;
		}

		std::random_device r;
		std::mt19937 generator(r());
		std::normal_distribution<float> distribution(0., 1.);

		for (long int j=0; j<nchans; j++)
		{
			if (snr2[j/fd] > threRFI2)
			{
				for (long int i=start[j/fd]; i<end[j/fd]; i++)
				{
					if (filltype == "mean")
						databuffer.buffer[i*nchans+j] = 0.;
					else
						databuffer.buffer[i*nchans+j] = distribution(generator);
				}
			}
		}

		//<0
		for (long int i=0; i<nsamples_ds; i++)
		{
			for (long int j=0; j<nchans_ds; j++)
			{
				buffer_ds[i*nchans_ds+j] = -buffer_ds[i*nchans_ds+j];
			}
		}

		PulsarX::kadane2D(boxsum.data(), start.data(), end.data(), buffer_ds.data(), nsamples_ds, nchans_ds);

		for (long int j=0; j<nchans_ds; j++)
		{
			wn[j] = end[j]-start[j]+1;

			float mean = 0.;
			float var = td*fd;

			if (wn[j] > wnlimit)
			{
				snr2[j] = boxsum[j]*boxsum[j]/(wn[j]*var);
			}
			else
			{
				snr2[j] = 0.;
			}

			start[j] *= td;
			end[j] += 1;
			end[j] *= td;
		}

		for (long int j=0; j<nchans; j++)
		{
			if (snr2[j/fd] > threRFI2)
			{
				for (long int i=start[j/fd]; i<end[j/fd]; i++)
				{
					if (filltype == "mean")
						databuffer.buffer[i*nchans+j] = 0.;
					else
						databuffer.buffer[i*nchans+j] = distribution(generator);
				}
			}
		}

		equalized = databuffer.equalized;

		databuffer.isbusy = true;

		BOOST_LOG_TRIVIAL(debug)<<"finished";

		return databuffer.get();
	}
	else
	{
		if (!databuffer.equalized)
		{
			BOOST_LOG_TRIVIAL(warning)<<"Warning: data is not equalize, kadaneF filter will not be performed"<<endl;
			return databuffer.get();
		}

		BOOST_LOG_TRIVIAL(debug)<<"perform kadane filter";

		if (closable) open();

		vector<float> bufferT(nchans*nsamples, 0.);

		transpose_pad<float>(&bufferT[0], &databuffer.buffer[0], nsamples, nchans);

		long int nsamples_ds = nsamples/td;
		long int nchans_ds = nchans/fd;

		vector<float> bufferT_ds(nchans_ds*nsamples_ds, 0.);

		//downsample
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(num_threads)
	#endif
		for (long int j=0; j<nchans_ds; j++)
		{
			for (long int k=0; k<fd; k++)
			{
				for (long int n=0; n<td; n++)
				{
					for (long int i=0; i<nsamples_ds; i++)       
					{
						bufferT_ds[j*nsamples_ds+i] += bufferT[(j*fd+k)*nsamples+(i*td+n)];
					}
				}
			}
		}

	#ifdef _OPENMP
		float *chdata_t = new float [num_threads*nsamples_ds];
		memset(chdata_t, 0, sizeof(float)*num_threads*nsamples_ds);

		std::vector<std::random_device> r(num_threads);
		std::vector<std::mt19937> generators;
		for (long int k=0; k<num_threads; k++) generators.emplace_back(std::mt19937(r[k]()));
		std::vector<std::normal_distribution<float>> distributions(num_threads, std::normal_distribution<float>(0., 1.));
	#else
		float *chdata_t = new float [nsamples_ds];
		memset(chdata_t, 0, sizeof(float)*nsamples_ds);

		std::random_device r;
		std::mt19937 generator(r());
		std::normal_distribution<float> distribution(0., 1.);
	#endif

		int wnlimit = widthlimit/tsamp/td;

	#ifdef _OPENMP
	#pragma omp parallel for num_threads(num_threads)
	#endif
		for (long int j=0; j<nchans_ds; j++)
		{
	#ifdef _OPENMP
			int thread_id = omp_get_thread_num();
			float *chdata = chdata_t+thread_id*nsamples_ds;
	#else
			float *chdata = chdata_t;
	#endif

			for (long int i=0; i<nsamples_ds; i++)
			{
				chdata[i] = bufferT_ds[j*nsamples_ds+i];
			}

			long int start, end;
			float snr2=0.;
			float boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
			long int wn = end-start+1;

			float mean = 0.;
			float var = td*fd;    
			if (wn > wnlimit)
			{
				snr2 = boxsum*boxsum/(wn*var);
			}
			else
			{
				snr2 = 0.;
			}

			start *= td;
			end += 1;
			end *= td;
			if (snr2 > threRFI2)
			{
				if (filltype == "mean")
				{
					for (long int k=0; k<fd; k++)
					{
						for (long int i=start; i<end; i++)
						{
							bufferT[(j*fd+k)*nsamples+i] = 0.;
						}
					}
				}
				else
				{
					for (long int k=0; k<fd; k++)
					{
						for (long int i=start; i<end; i++)
						{
	#ifdef _OPENMP
							bufferT[(j*fd+k)*nsamples+i] = distributions[thread_id](generators[thread_id]);
	#else
							bufferT[(j*fd+k)*nsamples+i] = distribution(generator);
	#endif
						}
					}
				}
			}

			//<0
			for (long int i=0; i<nsamples_ds; i++)
			{
				chdata[i] = -chdata[i];
			}

			snr2=0.;
			boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
			wn = end-start+1;

			if (wn > wnlimit)
			{
				snr2 = boxsum*boxsum/(wn*var);
			}
			else
			{
				snr2 = 0.;
			}

			start *= td;
			end += 1;
			end *= td;
			if (snr2 > threRFI2)
			{
				if (filltype == "mean")
				{
					for (long int k=0; k<fd; k++)
					{
						for (long int i=start; i<end; i++)
						{
							bufferT[(j*fd+k)*nsamples+i] = 0.;
						}
					}
				}
				else
				{
					for (long int k=0; k<fd; k++)
					{
						for (long int i=start; i<end; i++)
						{
	#ifdef _OPENMP
							bufferT[(j*fd+k)*nsamples+i] = distributions[thread_id](generators[thread_id]);
	#else
							bufferT[(j*fd+k)*nsamples+i] = distribution(generator);
	#endif
						}
					}
				}
			}
			
		}

		transpose_pad<float>(&buffer[0], &bufferT[0], nchans, nsamples);

		equalized = databuffer.equalized;

		databuffer.isbusy = false;
		isbusy = true;

		if (databuffer.closable) databuffer.close();

		delete [] chdata_t;

		BOOST_LOG_TRIVIAL(debug)<<"finished";

		return this;		
	}
#endif
}

DataBuffer<float> * RFI::kadaneT(DataBuffer<float> &databuffer, float threRFI2, double bandlimit, int td, int fd)
{
	if (!databuffer.equalized)
	{
		cerr<<"Error: data is not equalize"<<endl;
		return databuffer.get();
	}

	if (closable) open();

	long int nsamples_ds = nsamples/td;
	long int nchans_ds = nchans/fd;

	vector<float> buffer_ds(nsamples_ds*nchans_ds, 0.);

	//downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples_ds; i++)
	{
		for (long int n=0; n<td; n++)
		{
			for (long int k=0; k<fd; k++)
			{
				for (long int j=0; j<nchans_ds; j++)
				{
					buffer_ds[i*nchans_ds+j] += databuffer.buffer[(i*td+n)*nchans+j*fd+k];
				}
			}
		}
	}

	buffer = databuffer.buffer;

#ifdef _OPENMP
	float *tsdata_t = new float [num_threads*nchans_ds];
	memset(tsdata_t, 0, sizeof(float)*num_threads*nchans_ds);
#else
	float *tsdata_t = new float [nchans_ds];
	memset(tsdata_t, 0, sizeof(float)*nchans_ds);
#endif

	int chnlimit = abs(bandlimit/(frequencies[1]-frequencies[0])/fd);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples_ds; i++)
	{
#ifdef _OPENMP
		float *tsdata = tsdata_t+omp_get_thread_num()*nchans_ds;
#else
		float *tsdata = tsdata_t;
#endif

		for (long int j=0; j<nchans_ds; j++)
		{
			tsdata[j] = buffer_ds[i*nchans_ds+j];
		}

		long int start, end;
		float snr2=0.;
		float boxsum = kadane<float>(tsdata, nchans_ds, &start, &end);
		long int chn = end-start+1;

		float mean = 0.;
		float var = td*fd;    

		if (chn > chnlimit)
		{
			snr2 = boxsum*boxsum/(chn*var);
		}
		else
		{
			snr2 = 0.;
		}

		start *= fd;
		end += 1;
		end *= fd;
		if (chn > nchans_ds*0.8)
		{
			start = 0;
			end = nchans-1;
		}

		if (snr2 > threRFI2)
		{
			for (long int k=0; k<td; k++)
			{
				for (long int j=start; j<end; j++)
				{
					buffer[(i*td+k)*nchans+j] = 0.;
				}
			}
		}

		//<0
		for (long int j=0; j<nchans_ds; j++)
		{
			tsdata[j] = -tsdata[j];
		}

		snr2=0.;
		boxsum = kadane<float>(tsdata, nchans_ds, &start, &end);
		chn = end-start+1;

		if (chn > chnlimit)
		{
			snr2 = boxsum*boxsum/(chn*var);
		}
		else
		{
			snr2 = 0.;
		}

		start *= fd;
		end += 1;
		end *= fd;
		
		if (chn > nchans_ds*0.8)
		{
			start = 0;
			end = nchans-1;
		}

		if (snr2 > threRFI2)
		{
			for (long int k=0; k<td; k++)
			{
				for (long int j=start; j<end; j++)
				{
					buffer[(i*td+k)*nchans+j] = 0.;
				}
			}
		}
	}

	equalized = databuffer.equalized;

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	delete [] tsdata_t;

	return this;
}
