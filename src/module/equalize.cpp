/*
 * equalize.cpp
 *
 *  Created on: May 5, 2020
 *      Author: ypmen
 */

#include <string.h>
#include <cmath>

#include "equalize.h"
#include "dedisperse.h"
#include "logging.h"

using namespace std;

Equalize::Equalize(){}

Equalize::Equalize(const Equalize &equalize) : DataBuffer<float>(equalize)
{}

Equalize & Equalize::operator=(const Equalize &equalize)
{
    DataBuffer<float>::operator=(equalize);
    return *this;
}

Equalize::~Equalize(){}

void Equalize::prepare(DataBuffer<float> &databuffer)
{
    nsamples = databuffer.nsamples;
    nchans = databuffer.nchans;

    resize(nsamples, nchans);

    tsamp = databuffer.tsamp;
    frequencies = databuffer.frequencies;

    chmean.resize(nchans, 0.);
    chstd.resize(nchans, 0.);
}

DataBuffer<float> * Equalize::run(DataBuffer<float> &databuffer)
{
    if (databuffer.equalized)
    {
        return databuffer.get();
    }

    BOOST_LOG_TRIVIAL(debug)<<"perform noramlization";

    if (closable) open();

    fill(chmean.begin(), chmean.end(), 0.);
    fill(chstd.begin(), chstd.end(), 0.);

    for (long int i=0; i<nsamples; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            chmean[j] += databuffer.buffer[i*nchans+j];
            chstd[j] += databuffer.buffer[i*nchans+j]*databuffer.buffer[i*nchans+j];
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
            buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-chmean[j])/chstd[j];
        }
    }

    equalized = true;
    counter += nsamples;

    databuffer.isbusy = false;
    isbusy = true;

    if (databuffer.closable) databuffer.close();

    BOOST_LOG_TRIVIAL(debug)<<"finished";

    return this;
}
