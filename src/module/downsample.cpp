/*
 * downsample.cpp
 *
 *  Created on: May 5, 2020
 *      Author: ypmen
 */

#include <string.h>

#include "dedisperse.h"
#include "downsample.h"

using namespace std;

Downsample::Downsample()
{
    td = 1;
    fd = 1;
}

Downsample::Downsample(const Downsample &downsample) : DataBuffer<float>(downsample)
{
    td = downsample.td;
    fd = downsample.fd;
}

Downsample & Downsample::operator=(const Downsample &downsample)
{
    DataBuffer<float>::operator=(downsample);

    td = downsample.td;
    fd = downsample.fd;

    return *this;
}

Downsample::Downsample(int tds, int fds)
{
    td = tds;
    fd = fds;
}

Downsample::~Downsample(){}

void Downsample::prepare(DataBuffer<float> &databuffer)
{
    assert(databuffer.nsamples%td == 0);

    nsamples = databuffer.nsamples/td;
    nchans = databuffer.nchans/fd;
    resize(nsamples, nchans);

    tsamp = databuffer.tsamp*td;

    fill(frequencies.begin(), frequencies.end(), 0.);
    for (long int j=0; j<nchans; j++)
    {
        for (long int k=0; k<fd; k++)
        {
            frequencies[j] += databuffer.frequencies[j*fd+k];
        }
        frequencies[j] /= fd;
    }
}

void Downsample::run(DataBuffer<float> &databuffer)
{
    fill(buffer.begin(), buffer.end(), 0.);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {
        for (long int n=0; n<td; n++)
        {
            for (long int k=0; k<fd; k++)
            {
                for (long int j=0; j<nchans; j++)
                {
                    buffer[i*nchans+j] += databuffer.buffer[(i*td+n)*databuffer.nchans+j*fd+k];
                }
            }
        }
    }
    
    equalized = false;
    counter += nsamples;
}
