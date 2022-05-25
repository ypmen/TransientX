/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-05-13 11:59:40
 * @modify date 2022-05-13 11:59:40
 * @desc [description]
 */

#include "dedisperse.h"
#include "patch.h"
#include "utils.h"
#include "logging.h"

#ifdef __AVX2__
#include "avx2.h"
#endif

Patch::Patch()
{
	filltype = "none";
	width = 0.1;
	threshold = 5.;
	killrate = 0.;
}

Patch::~Patch(){}

void Patch::prepare(DataBuffer<float> &databuffer)
{
	nsamples = databuffer.nsamples;
	nchans = databuffer.nchans;

	resize(nsamples, nchans);

	tsamp = databuffer.tsamp;
	frequencies = databuffer.frequencies;
}

DataBuffer<float> * Patch::filter(DataBuffer<float> &databuffer)
{
	if (filltype != "mean" and filltype != "rand")
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform running median filter with time scale="<<width;

	std::vector<double> sstd(nsamples, 0);
	std::vector<bool> mask(nsamples, false);

	if (nchans % 8 == 0)
	{
		for (long int i=0; i<nsamples; i++)
		{
			float temp1 = 0.;
			float temp2 = 0.;
#ifndef __AVX2__
			for (long int j=0; j<nchans; j++)
			{
				temp1 += databuffer.buffer[i*nchans+j];
				temp2 += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
#else
			PulsarX::accumulate_mean_var3(temp1, temp2, databuffer.buffer.data() + i*nchans, nchans);
#endif
			temp1 /= nchans;
			temp2 /= nchans;
			temp2 -= temp1 * temp1;

			sstd[i] = temp2;
		}
	}
	else
	{
		for (long int i=0; i<nsamples; i++)
		{
			float temp1 = 0.;
			float temp2 = 0.;
			for (long int j=0; j<nchans; j++)
			{
				temp1 += databuffer.buffer[i*nchans+j];
				temp2 += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
			temp1 /= nchans;
			temp2 /= nchans;
			temp2 -= temp1 * temp1;

			sstd[i] = temp2;
		}
	}

	long int kill_count = 0;

	for (long int i=0; i<nsamples; i++)
	{
		if (sstd[i] == 0.)
		{
			mask[i] = true;
			if (i != 0)
				mask[i-1] = true;
			if (i != nsamples - 1)
				mask[i+1] = true;
			kill_count++;
		}
	}

	killrate = kill_count * 1. / nsamples;

#ifndef __AVX2__
	std::vector<double> chvar(nchans, 0.);
	std::vector<double> chmean(nchans, 0.);
#else
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> chvar(nchans, 0.);
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> chmean(nchans, 0.);
#endif
	long int cnt = 0;

	if (nchans % 4 == 0)
	{
		for (long int i=0; i<nsamples; i++)
		{
			if (mask[i]) continue;

#ifndef __AVX2__
			for (long int j=0; j<nchans; j++)
			{
				chmean[j] += databuffer.buffer[i*nchans+j];
				chvar[j] += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
#else
			PulsarX::accumulate_mean_var(chmean.data(), chvar.data(), databuffer.buffer.data()+i*nchans, nchans);
#endif
			cnt++;
		}
	}
	else
	{
		for (long int i=0; i<nsamples; i++)
		{
			if (mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				chmean[j] += databuffer.buffer[i*nchans+j];
				chvar[j] += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
			cnt++;
		}
	}

	for (long int j=0; j<nchans; j++)
	{
		chmean[j] /= cnt;
		chvar[j] /= cnt;
		chvar[j] -= chmean[j] * chmean[j];
	}

	if (filltype == "mean")
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples; i++)
		{
			if (!mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] = chmean[j];
			}
		}
	}
	else if (filltype == "rand")
	{
		std::vector<std::random_device> r(nchans);
		std::vector<std::mt19937> generators;
		std::vector<std::normal_distribution<float>> distributions;
		for (long int j=0; j<nchans; j++)
		{
			generators.emplace_back(std::mt19937(r[j]()));
			distributions.emplace_back(std::normal_distribution<float>(chmean[j], std::sqrt(chvar[j])));
		}

		for (long int i=0; i<nsamples; i++)
		{
			if (!mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] = distributions[j](generators[j]);
			}
		}
	}

	databuffer.equalized = false;
	counter + nsamples;

	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished"<<"("<<"killrate = "<<killrate<<")";

	return databuffer.get();
}

DataBuffer<float> * Patch::filter2(DataBuffer<float> &databuffer)
{
	if (filltype != "mean" and filltype != "rand")
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform running median filter with time scale="<<width;

	std::vector<double> szero(nsamples, 0);
	std::vector<double> s(nsamples, 0);
	std::vector<bool> mask(nsamples, false);

	if (nchans % 8 == 0)
	{
		for (long int i=0; i<nsamples; i++)
		{
#ifndef __AVX2__
			double temp = 0.;
			for (long int j=0; j<nchans; j++)
			{
				temp += databuffer.buffer[i*nchans+j];
			}
#else
			double temp = PulsarX::reduce(databuffer.buffer.data() + i*nchans, nchans);
#endif
			szero[i] = temp/nchans;
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
			szero[i] = temp/nchans;
		}
	}

	runMedian2(szero.data(), s.data(), nsamples, width/tsamp);

	std::vector<double> ssort(nsamples);
	for (long int i=0; i<nsamples; i++)
	{
		ssort[i] = (szero[i] - s[i])*(szero[i] - s[i]);
	}

	std::nth_element(ssort.begin(), ssort.begin()+ssort.size()/4, ssort.end(), std::less<float>());
	float q1 = ssort[ssort.size()/4];
	std::nth_element(ssort.begin(), ssort.begin()+ssort.size()/4, ssort.end(), std::greater<float>());
	float q3 =ssort[ssort.size()/4];
	float R = q3-q1;

	long int kill_count = 0;

	for (long int i=0; i<nsamples; i++)
	{
		if ((szero[i]-s[i]) < q1 - threshold * R ||
			(szero[i]-s[i]) > q3 + threshold * R)
		{
			mask[i] = true;
			kill_count++;
		}
	}

	killrate = kill_count * 1. / nsamples;

#ifndef __AVX2__
	std::vector<double> chvar(nchans, 0.);
	std::vector<double> chmean(nchans, 0.);
#else
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> chvar(nchans, 0.);
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> chmean(nchans, 0.);
#endif
	long int cnt = 0;

	if (nchans % 4 == 0)
	{
		for (long int i=0; i<nsamples; i++)
		{
			if (mask[i]) continue;

#ifndef __AVX2__
			for (long int j=0; j<nchans; j++)
			{
				chmean[j] += databuffer.buffer[i*nchans+j];
				chvar[j] += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
#else
			PulsarX::accumulate_mean_var(chmean.data(), chvar.data(), databuffer.buffer.data()+i*nchans, nchans);
#endif
			cnt++;
		}
	}
	else
	{
		for (long int i=0; i<nsamples; i++)
		{
			if (mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				chmean[j] += databuffer.buffer[i*nchans+j];
				chvar[j] += databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
			cnt++;
		}
	}

	for (long int j=0; j<nchans; j++)
	{
		chmean[j] /= cnt;
		chvar[j] /= cnt;
		chvar[j] -= chmean[j] * chmean[j];
	}

	if (filltype == "mean")
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples; i++)
		{
			if (!mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] = chmean[j];
			}
		}
	}
	else if (filltype == "rand")
	{
		std::vector<std::random_device> r(nchans);
		std::vector<std::mt19937> generators;
		std::vector<std::normal_distribution<float>> distributions;
		for (long int j=0; j<nchans; j++)
		{
			generators.emplace_back(std::mt19937(r[j]()));
			distributions.emplace_back(std::normal_distribution<float>(chmean[j], std::sqrt(chvar[j])));
		}

		for (long int i=0; i<nsamples; i++)
		{
			if (!mask[i]) continue;

			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] = distributions[j](generators[j]);
			}
		}
	}

	databuffer.equalized = false;
	counter + nsamples;

	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished"<<"("<<"killrate = "<<killrate<<")";

	return databuffer.get();
}