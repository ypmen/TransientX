/*
 * downsample.cpp
 *
 *  Created on: May 5, 2020
 *      Author: ypmen
 */

#include <string.h>

#include "dedisperse.h"
#include "downsample.h"
#include "logging.h"

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

	if (td != 1 or fd !=1)
	{
		std::vector<std::pair<std::string, std::string>> meta = {
			{"nsamples", std::to_string(databuffer.nsamples)},
			{"nchans", std::to_string(databuffer.nchans)},
			{"tsamp", std::to_string(databuffer.tsamp)},
			{"td", std::to_string(td)},
			{"fd", std::to_string(fd)}
		};
		format_logging("Downsampling Info", meta);
	}
}

DataBuffer<float> * Downsample::run(DataBuffer<float> &databuffer)
{
	if (td == 1 && fd ==1)
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform downsampling width td="<<td<<" fd="<<fd;

	if (closable)
		open();
	else
		std::fill(buffer.begin(), buffer.end(), 0.);

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

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return this;
}
