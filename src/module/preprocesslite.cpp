/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-06-01 18:08:09
 * @modify date 2021-06-01 18:08:09
 * @desc [description]
 */

#include <limits>
#include <random>
#include "preprocesslite.h"
#include "utils.h"
#include "dedisperse.h"
#include "logging.h"

#ifdef __AVX2__
#include "avx2.h"
#endif

void PreprocessLite::prepare(DataBuffer<float> &databuffer)
{
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

	closable = false;

	std::vector<std::pair<std::string, std::string>> meta = {
		{"nsamples", std::to_string(databuffer.nsamples)},
		{"nchans", std::to_string(databuffer.nchans)},
		{"tsamp", std::to_string(databuffer.tsamp)},
		{"td", std::to_string(td)},
		{"fd", std::to_string(fd)},
		{"IQR threshold", std::to_string(thresig)}
	};
	format_logging("Preprocess Info", meta);
}

DataBuffer<float> * PreprocessLite::run(DataBuffer<float> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform skewness-kurtosis filter with iqr threshold="<<thresig;

	if (closable) open();

	std::vector<float> chkurtosis(databuffer.nchans, 0.), chskewness(databuffer.nchans, 0.), chmean(databuffer.nchans, 0.), chstd(databuffer.nchans, 0.);

#ifndef __AVX2__
	std::vector<double> chmean1(databuffer.nchans, 0.), chmean2(databuffer.nchans, 0.), chmean3(databuffer.nchans, 0.), chmean4(databuffer.nchans, 0.), chcorr(databuffer.nchans, 0.), last_data(databuffer.nchans, 0.);
	for (long int i=0; i<databuffer.nsamples; i++)
	{
		for (long int j=0; j<databuffer.nchans; j++)
		{
			double tmp1 = databuffer.buffer[i*databuffer.nchans+j];
			double tmp2 = tmp1*tmp1;
			double tmp3 = tmp2*tmp1;
			double tmp4 = tmp2*tmp2;
			chmean1[j] += tmp1;
			chmean2[j] += tmp2;
			chmean3[j] += tmp3;
			chmean4[j] += tmp4;
			
			chcorr[j] += tmp1*last_data[j];
			last_data[j] = tmp1;
		}
	}
#else
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> chmean1(databuffer.nchans, 0.), chmean2(databuffer.nchans, 0.), chmean3(databuffer.nchans, 0.), chmean4(databuffer.nchans, 0.), chcorr(databuffer.nchans, 0.), last_data(databuffer.nchans, 0.);

	if (databuffer.nchans % 4 == 0)
	{
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			PulsarX::accumulate_mean1_mean2_mean3_mean4_corr1(chmean1.data(), chmean2.data(), chmean3.data(), chmean4.data(), chcorr.data(), last_data.data(), databuffer.buffer.data()+i*nchans, nchans);
		}
	}
	else
	{
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans; j++)
			{
				double tmp1 = databuffer.buffer[i*databuffer.nchans+j];
				double tmp2 = tmp1*tmp1;
				double tmp3 = tmp2*tmp1;
				double tmp4 = tmp2*tmp2;
				chmean1[j] += tmp1;
				chmean2[j] += tmp2;
				chmean3[j] += tmp3;
				chmean4[j] += tmp4;
				
				chcorr[j] += tmp1*last_data[j];
				last_data[j] = tmp1;
			}
		}
	}

#endif

	for (long int j=0; j<databuffer.nchans; j++)
	{
		chmean1[j] /= databuffer.nsamples;
		chmean2[j] /= databuffer.nsamples;
		chmean3[j] /= databuffer.nsamples;
		chmean4[j] /= databuffer.nsamples;

		chcorr[j] /= databuffer.nsamples-1;

		double tmp = chmean1[j]*chmean1[j];

		chmean[j] = chmean1[j];
		chstd[j] = chmean2[j]-tmp;
		
		if (chstd[j] > 0.)
		{
			chskewness[j] = chmean3[j]-3.*chmean2[j]*chmean1[j]+2.*tmp*chmean1[j];
			chkurtosis[j] = chmean4[j]-4.*chmean3[j]*chmean1[j]+6.*chmean2[j]*tmp-3.*tmp*tmp;

			chkurtosis[j] /= chstd[j]*chstd[j];
			chkurtosis[j] -= 3.;

			chskewness[j] /= chstd[j]*std::sqrt(chstd[j]);

			chcorr[j] -= tmp;
			chcorr[j] /= chstd[j];
		}
		else
		{
			chstd[j] = 1.;
			chkurtosis[j] = std::numeric_limits<float>::max();
			chskewness[j] = std::numeric_limits<float>::max();

			chcorr[j] = std::numeric_limits<float>::max();
		}

		chstd[j] = std::sqrt(chstd[j]);
	}

	/* calculate mean and std of chkurtosis and chskewness */
	std::vector<float> kurtosis_sort = chkurtosis;
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::less<float>());
	float kurtosis_q1 = kurtosis_sort[kurtosis_sort.size()/4];
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::greater<float>());
	float kurtosis_q3 =kurtosis_sort[kurtosis_sort.size()/4];
	float kurtosis_R = kurtosis_q3-kurtosis_q1;

	std::vector<float> skewness_sort = chskewness;
	std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::less<float>());
	float skewness_q1 = skewness_sort[skewness_sort.size()/4];
	std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::greater<float>());
	float skewness_q3 =skewness_sort[skewness_sort.size()/4];
	float skewness_R = skewness_q3-skewness_q1;

	std::vector<double> corr_sort(chcorr.begin(), chcorr.end());
	std::nth_element(corr_sort.begin(), corr_sort.begin()+corr_sort.size()/4, corr_sort.end(), std::less<float>());
	double corr_q1 = corr_sort[corr_sort.size()/4];
	std::nth_element(corr_sort.begin(), corr_sort.begin()+corr_sort.size()/4, corr_sort.end(), std::greater<float>());
	double corr_q3 = corr_sort[corr_sort.size()/4];
	double corr_R = corr_q3-corr_q1;

	long int kill_count = 0;

	std::vector<float> weights(databuffer.nchans, 0.);
	for (long int j=0; j<databuffer.nchans; j++)
	{
		if (chkurtosis[j]>=kurtosis_q1-thresig*kurtosis_R && \
			chkurtosis[j]<=kurtosis_q3+thresig*kurtosis_R && \
			chskewness[j]>=skewness_q1-thresig*skewness_R && \
			chskewness[j]<=skewness_q3+thresig*skewness_R && \
			chcorr[j]>=corr_q1-thresig*corr_R && \
			chcorr[j]<=corr_q3+thresig*corr_R)
		{
			weights[j] = 1.;
		}
		else
		{
			kill_count++;
		}
	}

	killrate = kill_count * 1. / databuffer.nchans;

	if (td == 1 && fd == 1)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int i=0; i<nsamples; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				buffer[i*nchans+j] = weights[j]*(databuffer.buffer[i*databuffer.nchans+j]-chmean[j])/chstd[j];
			}
		}
	}
	else if (fd == 1)
	{
		std::fill(buffer.begin(), buffer.end(), 0);
		for (long int i=0; i<nsamples; i++)
		{
			for (long int m=0; m<td; m++)
			{
				for (long int j=0; j<nchans; j++)
				{
					buffer[i*nchans+j] += weights[j]*(databuffer.buffer[(i*td+m)*databuffer.nchans+j]-chmean[j])/chstd[j];
				}
			}
		}
	}
	else
	{
		std::fill(buffer.begin(), buffer.end(), 0);
		for (long int i=0; i<nsamples; i++)
		{
			for (long int m=0; m<td; m++)
			{
				for (long int l=0; l<fd; l++)
				{
					for (long int j=0; j<nchans; j++)
					{
						buffer[i*nchans+j] += weights[j*fd+l]*(databuffer.buffer[(i*td+m)*databuffer.nchans+(j*fd+l)]-chmean[j*fd+l])/chstd[j*fd+l];
					}
				}
			}
		}
	}

	double stddev = std::sqrt(td*fd);
	if (filltype == "rand")
	{
#ifdef _OPENMP
		std::vector<std::random_device> r(num_threads);
		std::vector<std::mt19937> generators;
		for (long int k=0; k<num_threads; k++) generators.emplace_back(std::mt19937(r[k]()));
		std::vector<std::normal_distribution<float>> distributions(num_threads, std::normal_distribution<float>(0., stddev));

		#pragma omp parallel for num_threads(num_threads)
		for (long int i=0; i<nsamples; i++)
		{
			int thread_id = omp_get_thread_num();
			for (long int j=0; j<nchans; j++)
			{
				if (weights[j] == 0.)
					buffer[i*nchans+j] = distributions[thread_id](generators[thread_id]);
			}
		}
#else
		std::random_device r;
		std::mt19937 generator(r());
		std::normal_distribution<float> distribution(0., stddev);

		for (long int i=0; i<nsamples; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				if (weights[j] == 0.)
					buffer[i*nchans+j] = distribution(generator);
			}
		}
#endif
	}

	if (td == 1 && fd == 1)
		equalized = true;
	else
		equalized = false;

	counter += nsamples;

	databuffer.isbusy = false;
	isbusy = true;

	if (databuffer.closable) databuffer.close();

	BOOST_LOG_TRIVIAL(debug)<<"finished"<<"("<<"killrate = "<<killrate<<")";

	return this;
}