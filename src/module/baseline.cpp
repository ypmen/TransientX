/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-01-17 14:26:45
 * @modify date 2021-01-17 14:26:45
 * @desc [description]
 */

#ifdef __AVX2__
#include "avx2.h"
#endif

#include "baseline.h"
#include "dedisperse.h"
#include "utils.h"
#include "logging.h"

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

	if (int(width/tsamp) >= 3)
	{
		std::vector<std::pair<std::string, std::string>> meta = {
			{"nsamples", std::to_string(nsamples)},
			{"nchans", std::to_string(nchans)},
			{"tsamp", std::to_string(tsamp)},
			{"width", std::to_string(width)}
		};
		format_logging("Baseline Removal Info", meta);
	}
}

DataBuffer<float> * BaseLine::filter(DataBuffer<float> &databuffer)
{
	if (int(width/tsamp) < 3)
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform baseline removal with time scale="<<width;

#ifndef __AVX2__
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
			databuffer.buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];

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
			databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-chmean[j])/chstd[j];
		}
	}

#else
	if (nchans % 8 == 0)
	{
		vector<double, boost::alignment::aligned_allocator<double, 32>> xe(nchans, 0.);
		vector<double, boost::alignment::aligned_allocator<double, 32>> xs(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> alpha(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> beta(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> szero(nsamples, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> s(nsamples, 0.);
		double se = 0.;
		double ss = 0.;

		for (long int i=0; i<nsamples; i++)
		{
			double temp = PulsarX::reduce(databuffer.buffer.data()+i*nchans, nchans);
			szero[i] = temp/nchans;
		}

		runMedian2(szero.data(), s.data(), nsamples, width/tsamp);

		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::accumulate_mean(xe.data(), xs.data(), s[i], databuffer.buffer.data()+i*nchans, nchans);

			se += s[i];
			ss += s[i]*s[i];
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

		vector<double, boost::alignment::aligned_allocator<double, 32>> chmean(nchans, 0.);
		vector<double, boost::alignment::aligned_allocator<double, 32>> chstd(nchans, 0.);

		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::remove_baseline(databuffer.buffer.data()+i*nchans, databuffer.buffer.data()+i*nchans, alpha.data(), beta.data(), s[i], nchans);
			PulsarX::accumulate_mean_var(chmean.data(), chstd.data(), databuffer.buffer.data()+i*nchans, nchans);
		}

		vector<float, boost::alignment::aligned_allocator<float, 32>> chmeanf(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> chstdf_inv(nchans, 0.);
		
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

			chmeanf[j] = chmean[j];
			chstdf_inv[j] = 1./chstd[j];
		}

	#ifdef _OPENMP
	#pragma omp parallel for num_threads(num_threads)
	#endif
		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::normalize2(databuffer.buffer.data()+i*nchans, databuffer.buffer.data()+i*nchans, chmeanf.data(), chstdf_inv.data(), nchans);
		}
	}
	else
	{
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
				databuffer.buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];

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
				databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-chmean[j])/chstd[j];
			}
		}
	}
#endif

	databuffer.equalized = true;
	counter += nsamples;

	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return databuffer.get();
}

DataBuffer<float> * BaseLine::run(DataBuffer<float> &databuffer)
{
	if (int(width/tsamp) < 3)
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform baseline removal with time scale="<<width;

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

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return this;
}

DataBuffer<float> * BaseLine::filter2(DataBuffer<float> &databuffer)
{
	if (int(width/tsamp) < 3)
	{
		return databuffer.get();
	}

	BOOST_LOG_TRIVIAL(debug)<<"perform baseline removal with time scale="<<width;

#ifndef __AVX2__
	vector<double> xe(nchans, 0.);
	vector<double> xs(nchans, 0.);
	vector<double> alpha(nchans, 0.);
	vector<double> beta(nchans, 0.);
	vector<double> szero(nsamples, 0.);
	vector<double> sstdzero(nsamples, 0.);
	vector<double> s(nsamples, 0.);
	vector<double> sstd(nsamples, 0.);
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
			databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j]);
			sstdzero[i] +=  databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
		}
		sstdzero[i] /= nchans;
	}

	runMedian2(sstdzero.data(), sstd.data(), nsamples, width/tsamp);

	for (long int i=0; i<nsamples; i++)
	{
		double norm = sstd[i] == 0. ? 0. : 1./sstd[i];
		for (long int j=0; j<nchans; j++)
		{
			databuffer.buffer[i*nchans+j] *= norm;

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
			databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-chmean[j])/chstd[j];
		}
	}

#else
	if (nchans % 8 == 0)
	{
		vector<double, boost::alignment::aligned_allocator<double, 32>> xe(nchans, 0.);
		vector<double, boost::alignment::aligned_allocator<double, 32>> xs(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> alpha(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> beta(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> szero(nsamples, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> sstdzero(nsamples, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> s(nsamples, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> sstd(nsamples, 0.);
		double se = 0.;
		double ss = 0.;

		for (long int i=0; i<nsamples; i++)
		{
			double temp = PulsarX::reduce(databuffer.buffer.data()+i*nchans, nchans);
			szero[i] = temp/nchans;
		}

		runMedian2(szero.data(), s.data(), nsamples, width/tsamp);

		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::accumulate_mean(xe.data(), xs.data(), s[i], databuffer.buffer.data()+i*nchans, nchans);

			se += s[i];
			ss += s[i]*s[i];
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

		vector<double, boost::alignment::aligned_allocator<double, 32>> chmean(nchans, 0.);
		vector<double, boost::alignment::aligned_allocator<double, 32>> chstd(nchans, 0.);

		for (long int i=0; i<nsamples; i++)
		{
			sstdzero[i] = PulsarX::remove_baseline_reduce(databuffer.buffer.data()+i*nchans, databuffer.buffer.data()+i*nchans, alpha.data(), beta.data(), s[i], nchans);
			sstdzero[i] /= nchans;
		}

		runMedian2(sstdzero.data(), sstd.data(), nsamples, width/tsamp);

		for (long int i=0; i<nsamples; i++)
		{
			double norm = sstd[i] == 0. ? 0. : 1./sstd[i];
			PulsarX::accumulate_mean_var_scale(chmean.data(), chstd.data(), databuffer.buffer.data()+i*nchans, norm, nchans);
		}

		vector<float, boost::alignment::aligned_allocator<float, 32>> chmeanf(nchans, 0.);
		vector<float, boost::alignment::aligned_allocator<float, 32>> chstdf_inv(nchans, 0.);
		
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

			chmeanf[j] = chmean[j];
			chstdf_inv[j] = 1./chstd[j];
		}

	#ifdef _OPENMP
	#pragma omp parallel for num_threads(num_threads)
	#endif
		for (long int i=0; i<nsamples; i++)
		{
			PulsarX::normalize2(databuffer.buffer.data()+i*nchans, databuffer.buffer.data()+i*nchans, chmeanf.data(), chstdf_inv.data(), nchans);
		}
	}
	else
	{
		vector<double> xe(nchans, 0.);
		vector<double> xs(nchans, 0.);
		vector<double> alpha(nchans, 0.);
		vector<double> beta(nchans, 0.);
		vector<double> szero(nsamples, 0.);
		vector<double> sstdzero(nsamples, 0.);
		vector<double> s(nsamples, 0.);
		vector<double> sstd(nsamples, 0.);
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
				databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j]);
				sstdzero[i] +=  databuffer.buffer[i*nchans+j] * databuffer.buffer[i*nchans+j];
			}
			sstdzero[i] /= nchans;
		}

		runMedian2(sstdzero.data(), sstd.data(), nsamples, width/tsamp);

		for (long int i=0; i<nsamples; i++)
		{
			double norm = sstd[i] == 0. ? 0. : 1./sstd[i];
			for (long int j=0; j<nchans; j++)
			{
				databuffer.buffer[i*nchans+j] *= norm;

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
				databuffer.buffer[i*nchans+j] = (databuffer.buffer[i*nchans+j]-chmean[j])/chstd[j];
			}
		}
	}
#endif

	databuffer.equalized = true;
	counter += nsamples;

	databuffer.isbusy = true;

	BOOST_LOG_TRIVIAL(debug)<<"finished";

	return databuffer.get();
}