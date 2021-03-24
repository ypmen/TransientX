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

DataBuffer<float> * BaseLine::run(DataBuffer<float> &databuffer)
{
    if (int(width/tsamp) < 3)
    {
        return databuffer.get();
    }

    if (closable) open();

    vector<double> xe(nchans, 0.);
    vector<double> xs(nchans, 0.);
    vector<double> alpha(nchans, 0.);
    vector<double> beta(nchans, 0.);
    vector<double> szero(nsamples, 0.);
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

    double tmp = se*se-ss*nsamples;
    if (tmp != 0)
    {
        for (long int j=0; j<nchans; j++)
        {
            alpha[j] = (xe[j]*se-xs[j]*nsamples)/tmp;
            beta[j] = (xs[j]*se-xe[j]*ss)/tmp;
        }
    }

    vector<double> chmean(nchans, 0.);
    vector<double> chstd(nchans, 0.);

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

    databuffer.isbusy = false;
    isbusy = true;

    if (databuffer.closable) databuffer.close();

    return this;
}
