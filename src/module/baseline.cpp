/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-01-17 14:26:45
 * @modify date 2021-01-17 14:26:45
 * @desc [description]
 */

#include "baseline.h"
#include "dedisperse.h"
#include "utils.h"

BaseLine::BaseLine()
{
    width = 0.;
}

BaseLine::~BaseLine(){}

void BaseLine::prepare(DataBuffer<float> &databuffer)
{
    nsamples = databuffer.nsamples;
    nchans = databuffer.nchans;

    resize(nsamples, nchans);

    tsamp = databuffer.tsamp;
    frequencies = databuffer.frequencies;
}

void BaseLine::run(DataBuffer<float> &databuffer)
{
    if (int(width/tsamp) < 3)
    {
        buffer = databuffer.buffer;
        equalized = databuffer.equalized;
        counter += nsamples;
        return;
    }

    vector<float> xe(nchans, 0.);
    vector<float> xs(nchans, 0.);
    vector<float> alpha(nchans, 0.);
    vector<float> beta(nchans, 0.);
    vector<float> szero(nsamples, 0.);
    vector<float> s(nsamples, 0.);
    float se = 0.;
    float ss = 0.;

    for (long int i=0; i<nsamples; i++)
    {
        float temp = 0.;
        for (long int j=0; j<nchans; j++)
        {
            temp += databuffer.buffer[i*nchans+j];
        }
        szero[i] = temp/nchans;
    }

    runMedian2(szero.data(), s.data(), nsamples, width/tsamp);

    for (long int i=0; i<nsamples; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            xe[j] += databuffer.buffer[i*nchans+j];
        }

        se += s[i];
        ss += s[i]*s[i];
        
        for (long int j=0; j<nchans; j++)
        {
            xs[j] += databuffer.buffer[i*nchans+j]*s[i];
        }
    }

    float tmp = se*se-ss*nsamples;
    for (long int j=0; j<nchans; j++)
    {
        alpha[j] = (xe[j]*se-xs[j]*nsamples)/tmp;
        beta[j] = (xs[j]*se-xe[j]*ss)/tmp;
    }

    vector<float> chmean(nchans, 0.);
    vector<float> chstd(nchans, 0.);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];

            chmean[j] += buffer[i*nchans+j];
            chstd[j] += buffer[i*nchans+j]*buffer[i*nchans+j];
        }
    }

    for (long int j=0; j<nchans; j++)
    {
        chmean[j] /= nsamples;
        chstd[j] /= nsamples;
        chstd[j] -= chmean[j]*chmean[j];
        chstd[j] = sqrt(chstd[j]);
        if (chstd[j] == 0)
        {
            chstd[j] = 1;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = (buffer[i*nchans+j]-chmean[j])/chstd[j];
        }
    }

    equalized = true;
    counter += nsamples;
}
