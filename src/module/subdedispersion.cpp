/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-07-28 09:45:43
 * @modify date 2020-07-28 09:45:43
 * @desc [description]
 */

#include <iomanip>
#include <algorithm>
#include <assert.h>
#include "subdedispersion.h"

using namespace std;
using namespace RealTime;

Subband::Subband()
{
    inplace = false;

    counter = 0;
    ndump = 0;
    noverlap = 0;
    nchans = 0;
    nsub = 0;
    ndm_per_sub = 0;
    ndm = 0;
    nsamples = 0;
    tsamp = 0.;
}

Subband::~Subband(){}

void Subband::prepare()
{
    double fmax = *max_element(frequencies.begin(), frequencies.end());

    mxdelayn.resize(nsub*nchans*ndm_per_sub, 0);
    for (long int k=0; k<nsub; k++)
    {
        for (long int j=0; j<nchans; j++)
        {
            for (long int i=0; i<ndm_per_sub; i++)
            {
                mxdelayn[k*nchans*ndm_per_sub+j*ndm_per_sub+i] = round(SubbandDedispersion::dmdelay(vdm[k*ndm_per_sub+i], fmax, frequencies[j])/tsamp);
            }
        }
    }

    buffer.resize(nsub*nsamples*nchans, 0.);
    bufferT.resize(nsub*nchans*nsamples, 0.);
    buffertim.resize(nsub*ndm_per_sub*ndump, 0.);
    cachetim.resize(nsub*ndm_per_sub*noverlap, 0.);
    cachesub.resize(nsub*nchans*(noverlap+(nsamples-ndump)), 0.);
}

void Subband::run(vector<float> &data)
{
    long int nspace = nsamples-ndump;
    for (long int k=0; k<nsub; k++)
    {
        for (long int i=0; i<ndump; i++)
        {
            for (long int j=0; j<nchans; j++)
            {
                buffer[k*nsamples*nchans+(i+nspace)*nchans+j] = data[k*ndump*nchans+i*nchans+j];
            }
        }
    }

    for (long int k=0; k<nsub; k++)
    {
        transpose_pad<float>(&bufferT[0]+k*nchans*nsamples, &buffer[0]+k*nsamples*nchans, nsamples, nchans);
    }

    fill(buffertim.begin(), buffertim.end(), 0.);
    for (long int k=0; k<nsub; k++)
    {
        for (long int j=0; j<nchans; j++)
        {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
            for (long int l=0; l<ndm_per_sub; l++)
            {
                for (long int i=0; i<ndump; i++)
                {
                    buffertim[(k*ndm_per_sub+l)*ndump+i] += bufferT[(k*nchans+j)*nsamples+(i+mxdelayn[(k*nchans+j)*ndm_per_sub+l])];
                }
            }
        }
    }

    for (long int k=0; k<nsub; k++)
    {
        for (long int i=0; i<nspace; i++)
        {
            for (long int j=0; j<nchans; j++)
            {
                buffer[k*nsamples*nchans+i*nchans+j] = buffer[k*nsamples*nchans+(i+ndump)*nchans+j];
            }
        }
    }

    counter += ndump;
}

void Subband::cache()
{
    for (long int k=0; k<nsub; k++)
    {
        for (long int l=0; l<ndm_per_sub; l++)
        {
            for (long int i=0; i<noverlap; i++)
            {
                cachetim[(k*ndm_per_sub+l)*noverlap+i] = buffertim[(k*ndm_per_sub+l)*ndump+(i+ndump-noverlap)];
            }
        }
    }

    for (long int k=0; k<nsub; k++)
    {
        for (long int j=0; j<nchans; j++)
        {
            for (long int i=0; i<noverlap+nsamples-ndump; i++)
            {
                cachesub[(k*nchans+j)*(noverlap+nsamples-ndump)+i] = bufferT[(k*nchans+j)*nsamples+(i+ndump-noverlap)];
            }
        }
    }
}

void Subband::get_subdata(vector<float> &subdata, int idm, bool overlaped) const
{
    if (!overlaped)
    {
        subdata.resize(nchans*ndump, 0.);

        int isub = idm/ndm_per_sub;
        int isubdm = idm%ndm_per_sub;

        if (!inplace)
        {
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<ndump; i++)
                {
                    subdata[j*ndump+i] = bufferT[(isub*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
        }
        else
        {
            int isub_real = decodeisub[isub];
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<ndump; i++)
                {
                    subdata[j*ndump+i] = bufferT[(isub_real*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
        }
    }
    else
    {
        subdata.resize(nchans*(noverlap+ndump), 0.);

        int isub = idm/ndm_per_sub;
        int isubdm = idm%ndm_per_sub;

        if (!inplace)
        {
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<noverlap; i++)
                {
                    subdata[j*(noverlap+ndump)+i] = cachesub[(isub*nchans+j)*(noverlap+(nsamples-ndump))+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<ndump; i++)
                {
                    subdata[j*(noverlap+ndump)+(i+noverlap)] = bufferT[(isub*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
        }
        else
        {
            int isub_real = decodeisub[isub];
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<noverlap; i++)
                {
                    subdata[j*(noverlap+ndump)+i] = cachesub[(isub_real*nchans+j)*(noverlap+(nsamples-ndump))+i+mxdelayn[(isub_real*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
            for (long int j=0; j<nchans; j++)
            {
                for (long int i=0; i<ndump; i++)
                {
                    subdata[j*(noverlap+ndump)+(i+noverlap)] = bufferT[(isub_real*nchans+j)*nsamples+i+mxdelayn[(isub_real*nchans+j)*ndm_per_sub+isubdm]];
                }
            }
        }
    }
}

void Subband::get_timdata(vector<float> &timdata, int idm, bool overlaped) const
{
    if (!overlaped)
    {
        timdata.resize(ndump, 0.);
        if (!inplace)
        {
            for (long int i=0; i<ndump; i++)
            {
                timdata[i] = buffertim[idm*ndump+i];
            }
        }
        else
        {
            int idm_real = decodeidm[idm];
            for (long int i=0; i<ndump; i++)
            {
                timdata[i] = buffertim[idm_real*ndump+i];
            }
        }
    }
    else
    {
        timdata.resize(noverlap+ndump, 0.);
        if (!inplace)
        {
            for (long int i=0; i<noverlap; i++)
            {
                timdata[i] = cachetim[idm*noverlap+i];
            }
            for (long int i=0; i<ndump; i++)
            {
                timdata[i+noverlap] = buffertim[idm*ndump+i];
            }
        }
        else
        {
            int idm_real = decodeidm[idm];
            for (long int i=0; i<noverlap; i++)
            {
                timdata[i] = cachetim[idm_real*noverlap+i];
            }
            for (long int i=0; i<noverlap+ndump; i++)
            {
                timdata[i+noverlap] = buffertim[idm_real*ndump+i];
            }
        }
    }
}

/** Subband dedisperison */
SubbandDedispersion::SubbandDedispersion()
{
    mean = 0.;
    var = 0.;
    counter = 0;
    offset = 0;
    noverlap = 0;
    nsubband = 0;
    ndump = 0;
    dms = 0.;
    ddm = 0.;
    ndm = 0;
    overlap = 0.;
    nchans = 0;
    nsamples = 0;
    tsamp = 0.;

    nsub = 0;
}

SubbandDedispersion::~SubbandDedispersion(){}

void SubbandDedispersion::prepare(DataBuffer<float> &databuffer)
{
    nchans = databuffer.nchans;
    tsamp = databuffer.tsamp;
    frequencies = databuffer.frequencies;

    nsubband = round(sqrt(nchans));

    double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchans; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}

    fmap.resize(nchans, 0);
    fcnt.resize(nsubband, 0);
    frefsub.resize(nsubband, 0.);

    /** map the subband*/
    double df2 = abs(1./(fmin*fmin)-1./(fmax*fmax))/nsubband;
    double fref = frequencies[0];
    vector<double> freq(nsubband, 0.);
    long int m=0;
    for (long int i=0; i<nsubband; i++)
    {
    	double fc0 = 0;
    	fcnt[i] = 0;
    	while(m<nchans)
    	{
    		double fc = frequencies[m];
    		if (abs(1./(fc*fc)-1./(fref*fref))/(i+1)<=df2)
    		{
    			fmap[m] = i;
    			fc0 += fc;
    			fcnt[i]++;
                frefsub[i] = fc>frefsub[i] ? fc:frefsub[i];
    			m++;
    		}
    		else
    			break;
    	}
    	freq[i]=fc0/fcnt[i];
    }

    double maxsubdelayN = ceil(dmdelay(dms+ndm*ddm, fmax, fmin)/nsubband/tsamp);
    double maxdelayN = ceil(dmdelay(dms+ndm*ddm, fmax, fmin)/tsamp);

    nsamples = maxsubdelayN+ndump;
    //assert(nsamples>=2*(nsamples-ndump));

    buffer.resize(nsamples*nchans, 0.);
    bufferT.resize(nchans*nsamples, 0.);
    
    /** prepare the subband */
    double ddm_sub = ddm*nsubband;
    nsub = ceil((float)ndm/nsubband);
    mxdelayn.resize(nchans*nsub, 0);
    
    noverlap = overlap*ndump;

    sub.rootname = rootname;
    sub.ndump = ndump;
    sub.noverlap = noverlap;
    sub.nchans = nsubband;
    sub.nsub = nsub;
    sub.ndm_per_sub = nsubband;
    sub.ndm = nsub*nsubband;
    sub.nsamples = maxdelayN+ndump;
    //assert(sub.nsamples>=2*(sub.nsamples-sub.ndump));
    sub.tsamp = tsamp;

    sub.vdm.resize(sub.ndm, 0.);
    for (long int i=0; i<sub.ndm; i++)
    {
        sub.vdm[i] = dms + i*ddm;
    }
    sub.frequencies.resize(nsubband, 0.);
    sub.frequencies = frefsub;
    sub.fcnt.resize(nsubband, 0);
    sub.fcnt = fcnt;
    sub.prepare();

    /** calculate mxdelayn */
    for (long int j=0; j<nchans; j++)
    {
        for (long int k=0; k<nsub; k++)
        {
            double dm = dms + k*ddm_sub;
            mxdelayn[j*nsub+k] = round(dmdelay(dm, frefsub[fmap[j]], frequencies[j])/tsamp);
        }
    }

    buffersub.resize(nsubband*nsub*ndump, 0.);
    buffersubT.resize(nsubband*nsub*ndump, 0.);

    offset = (nsamples-ndump)+(sub.nsamples-sub.ndump)+sub.noverlap;
}

void SubbandDedispersion::run(DataBuffer<float> &databuffer, long int ns)
{
    assert(ns == ndump);

    int nspace = nsamples-ns;

    for (long int i=0; i<ns; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[(i+nspace)*nchans+j] = databuffer.buffer[i*nchans+j];
        }
    }

    transpose_pad<float>(&bufferT[0], &buffer[0], nsamples, nchans);

    fill(buffersub.begin(), buffersub.end(), 0);
    for (long int j=0; j<nchans; j++)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int k=0; k<nsub; k++)
        {
            for (long int i=0; i<ndump; i++)
            {
                buffersub[fmap[j]*nsub*ndump+k*ndump+i] += bufferT[j*nsamples+i+mxdelayn[j*nsub+k]];
            }
        }
    }

    transpose_pad<float>(&buffersubT[0], &buffersub[0], nsubband, nsub*ndump);

    sub.run(buffersubT);

    for (long int i=0; i<nspace; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = buffer[(i+ns)*nchans+j];
        }
    }

    counter += ndump;
}

void SubbandDedispersion::preparedump()
{
    double fmin = 1e6;
	double fmax = 0.;
    for (long int j=0; j<nchans; j++)
    {
        fmax = frequencies[j]>fmax? frequencies[j]:fmax;
        fmin = frequencies[j]<fmin? frequencies[j]:fmin;
    }

    for (long int k=0; k<ndm; k++)
    {
        double dm = sub.vdm[k];
        stringstream ss_dm;
        ss_dm << "DM" /*<< setw(8)*/ << setprecision(2) << fixed << setfill('0') << dm;
        string s_dm = ss_dm.str();

        ofstream outfile;
        outfile.open(rootname + "_" + s_dm + ".dat", ios::binary|ios::app);
        outfile.write((char *)(&dm), sizeof(double));
        outfile.write((char *)(&tsamp), sizeof(double));
        outfile.write((char *)(&fmin), sizeof(double));
        outfile.write((char *)(&fmax), sizeof(double));
        outfile.close();
    }
}

void SubbandDedispersion::rundump()
{
    for (long int k=0; k<ndm; k++)
    {
        double dm = sub.vdm[k];
        stringstream ss_dm;
        ss_dm << "DM" /*<< setw(8)*/ << setprecision(2) << fixed << setfill('0') << dm;
        string s_dm = ss_dm.str();

        ofstream outfile;
        outfile.open(rootname + "_" + s_dm + ".dat", ios::binary|ios::app);
        outfile.write((char *)(&sub.buffertim[0]+k*sub.ndump), sizeof(float)*sub.ndump);
        outfile.close();
    }
}
