/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-13 21:19:49
 * @modify date 2021-04-13 21:19:49
 * @desc [description]
 */

#include <fftw3.h>

#include <plotx.h>
namespace plt = PlotX;

#include "logging.h"
#include "dedisperse.h"
#include "candidate.h"
#include "archivewriter_tx.h"
#include "utils.h"
#include "hough_transform.h"

Candidate::Candidate()
{
	captured = false;
	mjd_start = 0.;

	id = 0;
	mjd = 0.;
	dm = 0.;
	width = 0.;
	width_orig = 0.;
	snr = 0.;
	fl = 0.;
	fh = 0.;
	planid = 0;

	dm_maxsnr = 0.;

	poltype = AABBCRCI;
	tbin = 0.;
	npol = 0;
	nchan = 0;
	nbin = 0;

	dms = 0.;
	ddm = 0.;
	ndm = 0;

	isdedispersed = false;
	isnormalized = false;
	maxdmid = 0;
	maxtid = 0;
}

Candidate::~Candidate(){}

void Candidate::dedisperse(bool coherent)
{
	if (coherent)
	{
		std::vector<double> delays(nchan, 0.);
		for (long int j=0; j<nchan; j++)
		{
			delays[j] = dmdelay(dm, *(std::max_element(frequencies.begin(), frequencies.end())), frequencies[j]);
		}

		std::vector<std::complex<float>> temp(npol*nchan*nbin);

#ifdef _OPENMP
		fftwf_init_threads();
		fftwf_plan_with_nthreads(num_threads);
 #endif
		fftwf_plan plan_r2c = fftwf_plan_dft_r2c_1d(nbin, data.data(), reinterpret_cast<fftwf_complex *>(temp.data()), FFTW_ESTIMATE);
		fftwf_plan plan_c2r = fftwf_plan_dft_c2r_1d(nbin, reinterpret_cast<fftwf_complex *>(temp.data()), data.data(), FFTW_ESTIMATE);

		float norm = 1. / nbin;

		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				fftwf_execute_dft_r2c(plan_r2c, data.data()+k*nchan*nbin+j*nbin, reinterpret_cast<fftwf_complex *>(temp.data()+k*nchan*nbin+j*nbin));
				for (long int i=0; i<nbin/2+1; i++)
				{
					temp[k*nchan*nbin+j*nbin+i] *= norm * std::complex<float>(std::cos(2.*M_PI/nbin*i*delays[j]/tbin), std::sin(2.*M_PI/nbin*i*delays[j]/tbin));
				}

				for (long int i=nbin/2+1; i<nbin; i++)
				{
					temp[k*nchan*nbin+j*nbin+i] = std::conj(temp[k*nchan*nbin+j*nbin+(nbin-i)]);
				}
			}
		}

		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				fftwf_execute_dft_c2r(plan_c2r, reinterpret_cast<fftwf_complex *>(temp.data()+k*nchan*nbin+j*nbin), data.data()+k*nchan*nbin+j*nbin);
			}
		}

		fftwf_destroy_plan(plan_r2c);
		fftwf_destroy_plan(plan_c2r);
#ifdef _OPENMP
		fftwf_cleanup_threads();
#endif
	}
	else
	{
		std::vector<int> delayn(nchan, 0.);
		for (long int j=0; j<nchan; j++)
		{
			delayn[j] = std::round(dmdelay(dm, *(std::max_element(frequencies.begin(), frequencies.end())), frequencies[j])/tbin);
		}

		for (long int k=0; k<npol; k++)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
			for (long int j=0; j<nchan; j++)
			{
				std::vector<float> temp(nbin, 0.);
				std::memcpy(temp.data(), data.data()+k*nchan*nbin+j*nbin, sizeof(float)*nbin);

				int delay = delayn[j];
				delay %= nbin;
				if (delay<0) delay += nbin;
				int nleft = nbin - delay;
				for (long int i=0; i<nleft; i++)
				{
					data[k*nchan*nbin+j*nbin+i] = temp[i+delay];
				}
				for (long int i=nleft; i<nbin; i++)
				{
					data[k*nchan*nbin+j*nbin+i] = temp[delay-nbin+i];
				}
			}
		}
	}

	isdedispersed = true;
}

void Candidate::dededisperse(bool coherent)
{
	if (coherent)
	{
		std::vector<double> delays(nchan, 0.);
		for (long int j=0; j<nchan; j++)
		{
			delays[j] = dmdelay(-dm, *(std::max_element(frequencies.begin(), frequencies.end())), frequencies[j]);
		}

		std::vector<std::complex<float>> temp(npol*nchan*nbin);

#ifdef _OPENMP
		fftwf_init_threads();
		fftwf_plan_with_nthreads(num_threads);
 #endif
		fftwf_plan plan_r2c = fftwf_plan_dft_r2c_1d(nbin, data.data(), reinterpret_cast<fftwf_complex *>(temp.data()), FFTW_ESTIMATE);
		fftwf_plan plan_c2r = fftwf_plan_dft_c2r_1d(nbin, reinterpret_cast<fftwf_complex *>(temp.data()), data.data(), FFTW_ESTIMATE);

		float norm = 1. / nbin;

		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				fftwf_execute_dft_r2c(plan_r2c, data.data()+k*nchan*nbin+j*nbin, reinterpret_cast<fftwf_complex *>(temp.data()+k*nchan*nbin+j*nbin));
				for (long int i=0; i<nbin/2+1; i++)
				{
					temp[k*nchan*nbin+j*nbin+i] *= norm * std::complex<float>(std::cos(2.*M_PI/nbin*i*delays[j]/tbin), std::sin(2.*M_PI/nbin*i*delays[j]/tbin));
				}

				for (long int i=nbin/2+1; i<nbin; i++)
				{
					temp[k*nchan*nbin+j*nbin+i] = std::conj(temp[k*nchan*nbin+j*nbin+(nbin-i)]);
				}
			}
		}

		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				fftwf_execute_dft_c2r(plan_c2r, reinterpret_cast<fftwf_complex *>(temp.data()+k*nchan*nbin+j*nbin), data.data()+k*nchan*nbin+j*nbin);
			}
		}

		fftwf_destroy_plan(plan_r2c);
		fftwf_destroy_plan(plan_c2r);
#ifdef _OPENMP
		fftwf_cleanup_threads();
#endif
	}
	else
	{
		std::vector<int> delayn(nchan, 0.);
		for (long int j=0; j<nchan; j++)
		{
			delayn[j] = std::round(dmdelay(-dm, *(std::max_element(frequencies.begin(), frequencies.end())), frequencies[j])/tbin);
		}

		for (long int k=0; k<npol; k++)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
			for (long int j=0; j<nchan; j++)
			{
				std::vector<float> temp(nbin, 0.);
				std::memcpy(temp.data(), data.data()+k*nchan*nbin+j*nbin, sizeof(float)*nbin);

				int delay = delayn[j];
				delay %= nbin;
				if (delay<0) delay += nbin;
				int nleft = nbin - delay;
				for (long int i=0; i<nleft; i++)
				{
					data[k*nchan*nbin+j*nbin+i] = temp[i+delay];
				}
				for (long int i=nleft; i<nbin; i++)
				{
					data[k*nchan*nbin+j*nbin+i] = temp[delay-nbin+i];
				}
			}
		}
	}

	isdedispersed = false;
}

void Candidate::shrink_to_fit(int nwidth, bool pow2bin, int factor)
{
	long int nbin_new = nwidth*width/tbin;

	if (pow2bin)
	{
		nbin_new = std::pow(2, std::ceil(std::log2(nbin_new)));
	}

	nbin_new = std::ceil(nbin_new/factor)*factor;

	if (nbin_new >= nbin) return;

	std::vector<float> data_new(npol*nchan*nbin_new, 0);
	for (long int k=0; k<npol; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin_new; i++)
			{
				data_new[k*nchan*nbin_new+j*nbin_new+i] = data[k*nchan*nbin+j*nbin+i];
			}
		}
	}

	std::swap(data, data_new);
	nbin = nbin_new;
}

void Candidate::save2ar(const std::string &template_file)
{
	std::string temp;

	TransientX::ArchiveWriter writer;
	writer.npol = npol;
	writer.nchan = nchan;
	writer.nbin = nbin;
	
	writer.template_file = template_file;
	
	temp = pngname;
	temp.erase(temp.find(".png"), 4);
	writer.rootname = "!"+temp;

	temp = s_beam;
	temp.erase(0, 4);
	writer.ibeam = std::stoi(temp);

	writer.mode = Integration::FOLD;

	writer.src_name = obsinfo["Source_name"];
	writer.ra = obsinfo["RA"];
	writer.dec = obsinfo["DEC"];

	writer.start_mjd = MJD(mjd_start);
	writer.dm = dm;
	writer.frequencies = frequencies;

	writer.prepare();

	writer.run(data, npol, nchan, nbin, nbin*tbin, nbin*tbin, (mjd-mjd_start)*86400.);

	writer.close();
}

void Candidate::sumif()
{
	if (npol == 1) return;

	std::vector<float> data_new(nchan*nbin, 0.);

	int nifs = std::min(npol, 2);
	for (long int k=0; k<nifs; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				data_new[j*nbin+i] += data[k*nchan*nbin+j*nbin+i];
			}
		}
	}

	std::swap(data_new, data);
	npol = 1;
}

void Candidate::AABBCRCI2IQUV()
{
	if (poltype == PolType::AABBCRCI)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				float I = data[0*nchan*nbin+j*nbin+i] + data[1*nchan*nbin+j*nbin+i];
				float Q = data[0*nchan*nbin+j*nbin+i] - data[1*nchan*nbin+j*nbin+i];
				float U = 2. * data[2*nchan*nbin+j*nbin+i];
				float V = 2. * data[3*nchan*nbin+j*nbin+i];

				data[0*nchan*nbin+j*nbin+i] = I;
				data[1*nchan*nbin+j*nbin+i] = Q;
				data[2*nchan*nbin+j*nbin+i] = U;
				data[3*nchan*nbin+j*nbin+i] = V;
			}
		}

		poltype = PolType::IQUV;
	}
}
	
void Candidate::IQUV2AABBCRCI()
{
	if (poltype == PolType::IQUV)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				float xx = 0.5 * (data[0*nchan*nbin+j*nbin+i] + data[1*nchan*nbin+j*nbin+i]);
				float yy = 0.5 * (data[0*nchan*nbin+j*nbin+i] - data[1*nchan*nbin+j*nbin+i]);
				float xy = 0.5 * data[2*nchan*nbin+j*nbin+i];
				float yx = 0.5 * data[3*nchan*nbin+j*nbin+i];
				
				data[0*nchan*nbin+j*nbin+i] = xx;
				data[1*nchan*nbin+j*nbin+i] = yy;
				data[2*nchan*nbin+j*nbin+i] = xy;
				data[3*nchan*nbin+j*nbin+i] = yx;
			}
		}

		poltype = PolType::AABBCRCI;
	}
}

void Candidate::downsample(int td, int fd)
{
	if (td == 1 && fd == 1) return;

	int nchan_new = nchan/fd;
	int nbin_new = nbin/td;

	std::vector<float> data_new(npol*nchan_new*nbin_new, 0.);
	
	for (long int k=0; k<npol; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan_new; j++)
		{
			for (long int jj=0; jj<fd; jj++)
			{
				for (long int ii=0; ii<td; ii++)
				{
					for (long int i=0; i<nbin_new; i++)
					{
						data_new[k*nchan_new*nbin_new+j*nbin_new+i] += data[k*nchan*nbin+(j*fd+jj)*nbin+(i*td+ii)];
					}
				}
			}
		}
	}

	std::vector<double> frequencies_new(nchan_new, 0.);
	for (long int j=0; j<nchan_new; j++)
	{
		for (long int jj=0; jj<fd; jj++)
		{
			frequencies_new[j] += frequencies[j*fd+jj];
		}
		frequencies_new[j] /= fd;
	}

	std::swap(data_new, data);
	std::swap(frequencies_new, frequencies);
	nchan = nchan_new;
	nbin = nbin_new;
	tbin *= td;

	isnormalized = false;
}

void Candidate::get_stats()
{
	mean.clear();
	var.clear();
	skewness.clear();
	kurtosis.clear();

	mean.resize(npol*nchan, 0.);
	var.resize(npol*nchan, 0.);
	skewness.resize(nchan, 0.);
	kurtosis.resize(nchan, 0.);
	autocorr1.resize(nchan, 0.);

	int nifs = std::min(npol, 2);

	std::vector<float> data_sumif(nchan*nbin, 0.);
	for (long int k=0; k<nifs; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				data_sumif[j*nbin+i] += data[k*nchan*nbin+j*nbin+i];
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchan; j++)
	{
		float boxsum = 0.;
		for (long int i=0; i<nbin/2; i++)
		{
			boxsum += data_sumif[j*nbin+i];
		}

		float min = boxsum;
		long int istart = 0;
		long int iend = nbin/2;

		for (long int i=0; i<nbin; i++)
		{
			boxsum -= data_sumif[j*nbin+i];
			boxsum += data_sumif[j*nbin+(i+nbin/2)%nbin];

			if (boxsum < min)
			{
				min = boxsum;
				istart = i+1;
				iend = nbin/2+i+1;
			}
		}

		for (long int k=0; k<npol; k++)
		{
			for (long int i=istart; i<iend; i++)
			{
				mean[k*nchan+j] += data[k*nchan*nbin+j*nbin+i%nbin];
				var[k*nchan+j] += data[k*nchan*nbin+j*nbin+i%nbin]*data[k*nchan*nbin+j*nbin+i%nbin];
			}
			mean[k*nchan+j] /= nbin/2;
			var[k*nchan+j] /= nbin/2;
			var[k*nchan+j] -= mean[k*nchan+j]*mean[k*nchan+j];
		}

		float tmp_mean = min/(nbin/2);
		float tmp_var = 0.;
		float last_data = 0;
		for (long int i=istart; i<iend; i++)
		{
			float tmp = data_sumif[j*nbin+i%nbin]-tmp_mean;
			float tmp2 = tmp*tmp;
			float tmp3 = tmp2*tmp;
			float tmp4 = tmp3*tmp;
			tmp_var += tmp2;
			skewness[j] += tmp3;
			kurtosis[j] += tmp4;
			autocorr1[j] += tmp*last_data;
			last_data = tmp;
		}

		tmp_var /= nbin/2;

		if (tmp_var <=0)
		{
			skewness[j] = 0.;
			kurtosis[j] = 0.;
			autocorr1[j] = 0.;
		}

		skewness[j] /= nbin/2;
		skewness[j] /= tmp_var*std::sqrt(tmp_var);

		kurtosis[j] /= nbin/2;
		kurtosis[j] /= tmp_var*tmp_var;
		kurtosis[j] -= 3.;

		autocorr1[j] /= nbin/2 -1;
		autocorr1[j] /= tmp_var;
	}
}

void Candidate::normalize()
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<npol*nchan; j++)
	{
		if (var[j] <= 0.)
		{
			for (long int i=0; i<nbin; i++)
			{
				data[j*nbin+i] -= mean[j];
			}
			mean[j] = 0.;
			var[j] = 0.;
		}
		else
		{
			float std = std::sqrt(var[j]);
			for (long int i=0; i<nbin; i++)
			{
				data[j*nbin+i] -= mean[j];
				data[j*nbin+i] /= std;
			}

			mean[j] = 0.;
			var[j] = 1.;
		}
	}

	isnormalized = true;
}

void Candidate::azap(float threshold)
{
	weights.resize(nchan, 1.);

	std::vector<float> kurtosis_sort = kurtosis;
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::less<float>());
	float kurtosis_q1 = kurtosis_sort[kurtosis_sort.size()/4];
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::greater<float>());
	float kurtosis_q3 = kurtosis_sort[kurtosis_sort.size()/4];
	float kurtosis_R = kurtosis_q3-kurtosis_q1;

	std::vector<float> skewness_sort = skewness;
	std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::less<float>());
	float skewness_q1 = skewness_sort[skewness_sort.size()/4];
	std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::greater<float>());
	float skewness_q3 = skewness_sort[skewness_sort.size()/4];
	float skewness_R = skewness_q3-skewness_q1;

	std::vector<float> autocorr1_sort = autocorr1;
	std::nth_element(autocorr1_sort.begin(), autocorr1_sort.begin()+autocorr1_sort.size()/4, autocorr1_sort.end(), std::less<float>());
	float autocorr1_q1 = autocorr1_sort[autocorr1_sort.size()/4];
	std::nth_element(autocorr1_sort.begin(), autocorr1_sort.begin()+autocorr1_sort.size()/4, autocorr1_sort.end(), std::greater<float>());
	float autocorr1_q3 = autocorr1_sort[autocorr1_sort.size()/4];
	float autocorr1_R = autocorr1_q3-autocorr1_q1;

	for (long int j=0; j<nchan; j++)
	{
		if ((kurtosis[j]<kurtosis_q1-threshold*kurtosis_R) ||
			(kurtosis[j]>kurtosis_q3+threshold*kurtosis_R) ||
			(skewness[j]<skewness_q1-threshold*skewness_R) ||
			(skewness[j]>skewness_q3+threshold*skewness_R) ||
			(autocorr1[j]<autocorr1_q1-threshold*autocorr1_R) ||
			(autocorr1[j]>autocorr1_q3+threshold*autocorr1_R))
		{
			weights[j] = 0.;
			var[j] = 0.;
			mean[j] = 0;
		}
	}

	for (long int k=0; k<npol; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				data[k*nchan*nbin+j*nbin+i] *= weights[j];
			}
		}
	}
}

void Candidate::zap(const std::vector<std::pair<double, double>> &zaplist)
{
	if (zaplist.empty()) return;

	weights.resize(nchan, 1.);

	for (long int j=0; j<nchan; j++)
	{
		for (auto k=zaplist.begin(); k!=zaplist.end(); ++k)
		{
			if (frequencies[j]>=(*k).first and frequencies[j]<=(*k).second)
			{
				weights[j] = 0;
			}
		}
	}

	for (long int k=0; k<npol; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				data[k*nchan*nbin+j*nbin+i] *= weights[j];
			}
		}
	}
}

void Candidate::clip(int td, int fd, float threshold)
{
	if (!isnormalized)
	{
		BOOST_LOG_TRIVIAL(warning)<<"data is not normalized, so clip filter is not performed";
		return;
	}

	long int nsamples_ds = nbin/td;
	long int nchans_ds = nchan/fd;

	vector<float> data_ds(nchans_ds*nsamples_ds, 0.);

	// downsample
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
					data_ds[j*nsamples_ds+i] += data[(j*fd+k)*nbin+(i*td+n)];
				}
			}
		}
	}

	float tempthre = threshold*std::sqrt(td*fd);

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
					if (data_ds[j*nsamples_ds+i] > tempthre)
					{
						data[(j*fd+k)*nbin+(i*td+n)] = 0.;
					}
				}
			}
		}
	}
}

void Candidate::kadaneF(int td, int fd, float threshold, int nwidth)
{
	if (!isnormalized)
	{
		BOOST_LOG_TRIVIAL(warning)<<"data is not normalized, so kadane filter is not performed";
		return;
	}

	long int nsamples_ds = nbin/td;
	long int nchans_ds = nchan/fd;

	vector<float> data_ds(nchans_ds*nsamples_ds, 0.);

	// downsample
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
					data_ds[j*nsamples_ds+i] += data[(j*fd+k)*nbin+(i*td+n)];
				}
			}
		}
	}

	// kadane filter
#ifdef _OPENMP
	std::vector<float> chdata_t(num_threads*nsamples_ds, 0.);
#else
	std::vector<float> chdata_t(nsamples_ds, 0.);
#endif

	int wnlimit = nwidth*width/tbin/td;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchans_ds; j++)
	{
#ifdef _OPENMP
		float *chdata = chdata_t.data()+omp_get_thread_num()*nsamples_ds;
#else
		float *chdata = chdata_t.data();
#endif

		for (long int i=0; i<nsamples_ds; i++)
		{
			chdata[i] = data_ds[j*nsamples_ds+i];
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
		if (snr2 > threshold*threshold)
		{
			for (long int k=0; k<fd; k++)
			{
				for (long int i=start; i<end; i++)
				{
					data[(j*fd+k)*nbin+i] = 0.;
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
		if (snr2 > threshold*threshold)
		{
			for (long int k=0; k<fd; k++)
			{
				for (long int i=start; i<end; i++)
				{
					data[(j*fd+k)*nbin+i] = 0.;
				}
			}
		}
	}
}

void Candidate::zdot()
{
	if (npol != 1) return;

	std::vector<double> s(nbin, 0.);
	// normalize
	for (long int j=0; j<nchan; j++)
	{
		double tmpmean=0., tmpstd=0.;
		for (long int i=0; i<nbin; i++)
		{
			tmpmean += data[j*nbin+i];
			tmpstd += data[j*nbin+i]*data[j*nbin+i];
		}

		tmpmean /= nbin;
		tmpstd /= nbin;
		tmpstd -= tmpmean*tmpmean;
		tmpstd = std::sqrt(tmpstd);
		if (tmpstd == 0.) tmpstd = 1.;

		for (long int i=0; i<nbin; i++)
		{
			data[j*nbin+i] -= tmpmean;
			data[j*nbin+i] /= tmpstd;

			s[i] += data[j*nbin+i];
		}
	}

	for (long int i=0; i<nbin; i++)
	{
		s[i] /= nchan;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchan; j++)
	{
		double xe = 0.;
		double xs = 0.;
		double se = 0.;
		double ss = 0.;
		double alpha = 0.;
		double beta = 0.;
		for (long int i=0; i<nbin; i++)
		{
			xe += data[j*nbin+i];
			xs += data[j*nbin+i]*s[i];
			se += s[i];
			ss += s[i]*s[i];
		}

		double tmp = se*se-ss*nbin;
		if (tmp != 0)
		{
			alpha = (xe*se-xs*nbin)/tmp;
			beta = (xs*se-xe*ss)/tmp;
		}

		for (long int i=0; i<nbin; i++)
		{
			data[j*nbin+i] -= alpha*s[i]+beta;
		}
	}
}

void Candidate::optimize()
{
	/* calculate dm range */
	double tmp = dmdelay(1., *std::max_element(frequencies.begin(), frequencies.end()), *std::min_element(frequencies.begin(), frequencies.end()));
	double dmrange = (nbin*tbin)/tmp;
	dms = std::max(0., dm-dmrange/2.);
	//dms = dm-dmrange/2.;
	double dme = dm+dmrange/2.;

	int nchans = std::pow(2, std::ceil(std::log2(nchan)));
	int nsamples = nbin;

	ndm = nchans;
	ddm = (dme-dms)/ndm;

	int nifs = std::min(npol, 2);

	std::vector<float> buffer(nchans*nsamples, 0.);
	std::vector<double> frequenciesP(nchans, 0.);
	double foff = frequencies[1]-frequencies[0];
	if (foff < 0.)
	{
		for (long int k=0; k<nifs; k++)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
			for (long int j=0; j<nchan; j++)
			{
				for (long int i=0; i<nbin; i++)
				{
					buffer[j*nbin+i] += data[k*nchan*nbin+j*nbin+i];
				}
			}
		}

		for (long int j=0; j<nchan; j++)
		{
			frequenciesP[j] = frequencies[j];
		}
		for (long int j=nchan; j<nchans; j++)
		{
			frequenciesP[j] = frequencies[0]+j*foff;
		}
	}
	else
	{
		for (long int k=0; k<nifs; k++)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
			for (long int j=0; j<nchan; j++)
			{
				for (long int i=0; i<nbin; i++)
				{
					buffer[j*nbin+i] += data[k*nchan*nbin+(nchan-1-j)*nbin+i];
				}
			}
		}

		for (long int j=0; j<nchan; j++)
		{
			frequenciesP[j] = frequencies[nchan-1-j];
		}
		for (long int j=nchan; j<nchans; j++)
		{
			frequenciesP[j] = frequencies[0]+(nchan-1-j)*foff;
		}
	}

	std::function<int(int, int, int)> get_shift;
	if (isdedispersed == false)
		get_shift = [&](int x0, int x1, int p) {return std::round(dmdelay(dms+p*ddm, frequenciesP[x0], frequenciesP[x1])/tbin);};
	else
		get_shift = [&](int x0, int x1, int p) {return std::round(dmdelay(dms+p*ddm-dm, frequenciesP[x0], frequenciesP[x1])/tbin);};

	hough_transform(buffer, dmtmap, nchans, nsamples, get_shift);
}

void Candidate::matched_filter(float snrloss)
{
	// calculate boxcar width series
	float wfactor = 1./((1.-snrloss)*(1.-snrloss));
	std::vector<float> vwn;
	vwn.push_back(1);
	while (true)
	{
		int tmp_wn1 = vwn.back();
		int tmp_wn2 = tmp_wn1*wfactor;
		tmp_wn2 = tmp_wn2<tmp_wn1+1 ? tmp_wn1+1:tmp_wn2;
		if (tmp_wn2 > nbin/2) break;
		vwn.push_back(tmp_wn2);
	}
	int nbox = vwn.size();

	// boxcar filter
	dmtsnrmap.resize(ndm*nbin, 0.);
	dmtwidmap.resize(ndm*nbin, 0.);
	profile.resize(nbin, 0.);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<ndm; j++)
	{
		std::vector<float> tim(dmtmap.begin()+j*nbin, dmtmap.begin()+(j+1)*nbin);

		double tmpmean = 0.;
		double tmpvar = 0.;

		get_mean_var<std::vector<float>::iterator>(tim.begin(), tim.size(), tmpmean, tmpvar);

		std::vector<float> csump(2*nbin, 0.);
		float csum = 0.;
		for (long int i=0; i<nbin; i++)
		{
			csum += tim[i] - tmpmean;
			csump[i] = csum;
		}
		for (long int i=nbin; i<2*nbin; i++)
		{
			csum += tim[i-nbin] - tmpmean;
			csump[i] = csum;
		}

		std::vector<float> vS(nbin, 0.);
		std::vector<int> vwn_maxS(nbin, 0.);

		if (tmpvar > 0)
		{
			for (long int k=0; k<nbox; k++)
			{
				int wn = vwn[k];
				int wl = wn/2+1;
				int wh = (wn-1)/2;

				float temp = sqrt(1./(wn*tmpvar));

				for (long int i=wl; i<nbin; i++)
				{
					float boxsum = csump[i+wh]-csump[i-wl];
					float S = boxsum*temp;
					if (S>=vS[i])
					{
						vS[i] = S;
						vwn_maxS[i] = wn;
					}
				}
				for (long int i=nbin; i<nbin+wl; i++)
				{
					float boxsum = csump[i+wh]-csump[i-wl];
					float S = boxsum*temp;
					if (S>=vS[i-nbin])
					{
						vS[i-nbin] = S;
						vwn_maxS[i-nbin] = wn;
					}
				}
			}
		}

		std::copy(vS.begin(), vS.end(), dmtsnrmap.begin()+j*nbin);
		std::copy(vwn_maxS.begin(), vwn_maxS.end(), dmtwidmap.begin()+j*nbin);
	}

	maxdmid = 0;
	maxtid = 0;
	snr = dmtsnrmap[0];
	for (long int j=0; j<ndm; j++)
	{
		for (long int i=0; i<nbin; i++)
		{
			if (dmtsnrmap[j*nbin+i] > snr)
			{
				snr = dmtsnrmap[j*nbin+i];
				maxdmid = j;
				maxtid = i;
			}
		}
	}

	snr = dmtsnrmap[maxdmid*nbin+maxtid];
	width = dmtwidmap[maxdmid*nbin+maxtid]*tbin;
	dm_maxsnr = dms+maxdmid*ddm;
	mjd = mjd_start+(maxtid*tbin)/86400.;

	std::copy(dmtmap.begin()+maxdmid*nbin, dmtmap.begin()+(maxdmid+1)*nbin, profile.begin());
}

void Candidate::peak_search(float snrloss)
{
	// calculate boxcar width series
	float wfactor = 1./((1.-snrloss)*(1.-snrloss));
	std::vector<float> vwn;
	vwn.push_back(1);
	while (true)
	{
		int tmp_wn1 = vwn.back();
		int tmp_wn2 = tmp_wn1*wfactor;
		tmp_wn2 = tmp_wn2<tmp_wn1+1 ? tmp_wn1+1:tmp_wn2;
		if (tmp_wn2 > nbin/2) break;
		vwn.push_back(tmp_wn2);
	}
	int nbox = vwn.size();

	// boxcar filter
	dmtsnrmap.resize(nbin, 0.);
	dmtwidmap.resize(nbin, 0.);
	profile.resize(nbin, 0.);

	if (!isdedispersed)
		dedisperse();

	int nifs = std::min(npol, 2);
	for (long int k=0; k<nifs; k++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				profile[i] += data[k*nchan*nbin+j*nbin+i];
			}
		}
	}

	std::vector<float> tim(profile.begin(), profile.end());

	double tmpmean = 0.;
	double tmpvar = 0.;

	get_mean_var<std::vector<float>::iterator>(tim.begin(), tim.size(), tmpmean, tmpvar);

	std::vector<float> csump(2*nbin, 0.);
	float csum = 0.;
	for (long int i=0; i<nbin; i++)
	{
		csum += tim[i] - tmpmean;
		csump[i] = csum;
	}
	for (long int i=nbin; i<2*nbin; i++)
	{
		csum += tim[i-nbin] - tmpmean;
		csump[i] = csum;
	}

	std::vector<float> vS(nbin, 0.);
	std::vector<int> vwn_maxS(nbin, 0.);

	if (tmpvar > 0)
	{
		for (long int k=0; k<nbox; k++)
		{
			int wn = vwn[k];
			int wl = wn/2+1;
			int wh = (wn-1)/2;

			float temp = sqrt(1./(wn*tmpvar));

			for (long int i=wl; i<nbin; i++)
			{
				float boxsum = csump[i+wh]-csump[i-wl];
				float S = boxsum*temp;
				if (S>=vS[i])
				{
					vS[i] = S;
					vwn_maxS[i] = wn;
				}
			}
			for (long int i=nbin; i<nbin+wl; i++)
			{
				float boxsum = csump[i+wh]-csump[i-wl];
				float S = boxsum*temp;
				if (S>=vS[i-nbin])
				{
					vS[i-nbin] = S;
					vwn_maxS[i-nbin] = wn;
				}
			}
		}
	}

	std::copy(vS.begin(), vS.end(), dmtsnrmap.begin());
	std::copy(vwn_maxS.begin(), vwn_maxS.end(), dmtwidmap.begin());

	maxtid = 0;
	snr = dmtsnrmap[0];
	for (long int i=0; i<nbin; i++)
	{
		if (dmtsnrmap[i] > snr)
		{
			snr = dmtsnrmap[i];
			maxtid = i;
		}
	}

	snr = dmtsnrmap[maxtid];
	width = dmtwidmap[maxtid]*tbin;
	dm_maxsnr = dm;
	mjd = mjd_start+(maxtid*tbin)/86400.;
}

void Candidate::save2png(const std::string &rootname, float threS)
{
	std::string basename = pngname;
	basename.erase(basename.find(".png"), 4);

	stringstream ss_toa;
	ss_toa << setprecision(15) << fixed << mjd;
	string s_toa = ss_toa.str();

	stringstream ss_date;
	ss_date << setprecision(10) << fixed << std::stold(obsinfo["Date"]);
	string s_date = ss_date.str();

	string s_ibeam = obsinfo["Beam"];

	std::stringstream ss_snr;
	ss_snr<<fixed<<setprecision(1)<<snr;
	std::string s_snr = ss_snr.str();

	std::stringstream ss_width;
	ss_width<<fixed<<setprecision(1)<<width*1000;
	std::string s_width = ss_width.str();

	std::stringstream ss_dm;
	ss_dm<<fixed<<setprecision(1)<<dm_maxsnr;
	std::string s_dm = ss_dm.str();

	std::string s_ymw16_maxdm;
	std::stringstream ss_ymw16_maxdm;
	ss_ymw16_maxdm<<fixed<<setprecision(1)<<stod(obsinfo["MaxDM_YMW16"]);
	s_ymw16_maxdm = ss_ymw16_maxdm.str();

	double ymw16_dist = get_dist_ymw16(stod(obsinfo["GL"]), stod(obsinfo["GB"]), dm);
	std::stringstream ss_dist;
	ss_dist<<fixed<<setprecision(1)<<ymw16_dist;
	std::string s_dist = ss_dist.str();

	stringstream ss_id;
	ss_id << setw(2) << setfill('0') << planid;
	string s_id = ss_id.str();

	stringstream ss_k;
	ss_k << setw(2) << setfill('0') << id;
	string s_k = ss_k.str();

	string src_name = obsinfo["Source_name"];

	string figname = basename+"_replot.png";

	long int nsamples = nbin;
	long int nchans = nchan;

	int wn = width/tbin;

	std::vector<float> mxft(nchans*nsamples, 0.);
	for (long int k=0; k<npol; k++)
	{
		if (k<2)
		{
			for (long int j=0; j<nchan; j++)
			{
				for (long int i=0; i<nbin; i++)
				{
					mxft[j*nsamples+i] += data[k*nchan*nbin+j*nbin+i];
				}
			}
		}
	}

	std::vector<float> vt, vf, vpf, vpfcum;
	vt.resize(nsamples, 0.);
	vf.resize(nchans, 0.);
	vpf.resize(nchans, 0.);
	vpfcum.resize(nchans, 0.);

	for (long int i=0; i<nsamples; i++)
	{
		vt[i] = i*tbin;
	}

	vector<float> vfnorm(nchans, 0.);

	float fmax = 0.;
	float fmin = 1e6;
	for (long int j=0; j<nchans; j++)
	{
		vf[j] = frequencies[j];
		vfnorm[j] = frequencies[j];
		vpf[j] = 0.;
		vpfcum[j] = 0.;
		fmax = vf[j]>fmax ? vf[j]:fmax;
		fmin = vf[j]<fmin ? vf[j]:fmin;
	}

	for (long int j=0; j<nchans; j++)
	{
		for (long int i=0; i<nsamples; i++)
		{
			if (i>maxtid-width/tbin and i<maxtid+width/tbin)
			{
				vpf[j] += mxft[j*nsamples+i];
			}
		}
	}

	float cum = 0.;
	if (vf[1] > vf[0])
	{
		for (long int j=0; j<nchans; j++)
		{
			cum += vpf[j];
			vpfcum[j] = cum;
		}
	}
	else
	{
		for (long int j=nchans-1; j>=0; j--)
		{
			cum += vpf[j];
			vpfcum[j] = cum;
		}
	}

	long int tds = wn/8;
	tds = tds<1 ? 1:tds;
	tds = tds>nsamples ? nsamples:tds;
	std::vector<float> vt_down(nsamples/tds, 0.);
	long int fds = nchans/64;
	fds = fds<1 ? 1:fds;
	fds = fds>nchans ? nchans:fds;
	std::vector<float> vf_down(nchans/fds, 0.);
	std::vector<float> vfnorm_down(nchans/fds, 0.);
	std::vector<float> vpfcum_down(nchans/fds, 0.);
	std::vector<float> mxft_down((nchans/fds)*(nsamples/tds), 0.);
	for (long int i=0; i<nsamples/tds; i++)
	{
		for (long int ii=0; ii<tds; ii++)
		{
			vt_down[i] += vt[i*tds+ii];
		}
		vt_down[i] /= tds;
	}
	for (long int j=0; j<nchans/fds; j++)
	{
		for (long int jj=0; jj<fds; jj++)
		{
			vf_down[j] += vf[j*fds+jj];
			vfnorm_down[j] += vfnorm[j*fds+jj];
			vpfcum_down[j] += vpfcum[j*fds+jj];
		}
		vf_down[j] /= fds;
		vfnorm_down[j] /= fds;
		vpfcum_down[j] /= fds;
	}

	for (long int j=0; j<nchans/fds; j++)
	{
		for (long int jj=0; jj<fds; jj++)
		{
			for (long int i=0; i<nsamples/tds; i++)
			{
				for (long int ii=0; ii<tds; ii++)
				{
					mxft_down[j*(nsamples/tds)+i] += mxft[(j*fds+jj)*nsamples+(i*tds+ii)];
				}
			}
		}
	}

	std::vector<float> vp(profile.begin(), profile.end());

	float fl = 0.;
	float fh = 0.;

	long int fstart=0;
	long int fend=0;
	kadane<float>(vpf.data(), vpf.size(), &fstart, &fend);
	fl = vfnorm[fstart];
	fh = vfnorm[fend];

	float vp_down_std = 0.;
	float vp_down_mean = 0.;
	std::vector<float> vp_down(nsamples/tds, 0.);
	for (long int i=0; i<nsamples/tds; i++)
	{
		for (long int ii=0; ii<tds; ii++)
		{
			vp_down[i] += vp[i*tds+ii];
		}
		vp_down[i] /= tds;

		vp_down_mean += vp_down[i];
		vp_down_std += vp_down[i]*vp_down[i];
	}
	vp_down_mean /= nsamples/tds;
	vp_down_std /= nsamples/tds;
	vp_down_std -= vp_down_mean*vp_down_mean;
	vp_down_std = std::sqrt(vp_down_std);
	for (long int i=0; i<nsamples/tds; i++)
	{
		vp_down[i] -= vp_down_mean;
		vp_down[i] /= vp_down_std;
	}

	std::vector<float> mxdmt_down;
	std::vector<float> vdm_down;
	std::vector<float> vSdm_down;
	std::vector<float> vdm;

	if (ndm > 1)
	{
		std::vector<float> mxdmt(ndm*nsamples, 0.);
		vdm.resize(ndm, 0.);
		std::vector<float> vSdm(ndm, 0.);
	
		for (long int j=0; j<ndm; j++)
		{
			vdm[j] = dms+j*ddm;
			vSdm[j] = 0.;
		}
		for (long int j=0; j<ndm; j++)
		{
			for (long int i=0; i<nsamples; i++)
			{
				mxdmt[j*nsamples+i] = dmtsnrmap[j*nsamples+i];
				vSdm[j] = mxdmt[j*nsamples+i]>vSdm[j] ? mxdmt[j*nsamples+i]:vSdm[j];
			}
		}

		int dmds = ndm/64;
		dmds = dmds<1 ? 1:dmds;
		dmds = dmds>ndm ? ndm:dmds;

		mxdmt_down.resize((ndm/dmds)*(nsamples/tds), 0.);
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int jj=0; jj<dmds; jj++)
			{
				for (long int i=0; i<nsamples/tds; i++)
				{
					for (long int ii=0; ii<tds; ii++)
					{
						mxdmt_down[j*(nsamples/tds)+i] += mxdmt[(j*dmds+jj)*nsamples+i*tds+ii];
					}
				}
			}
		}
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int i=0; i<nsamples/tds; i++)
			{
				mxdmt_down[j*(nsamples/tds)+i] /= tds*dmds;
			}
		}

		vdm_down.resize(ndm/dmds, 0.);
		vSdm_down.resize(ndm/dmds, 0.);
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int jj=0; jj<dmds; jj++)
			{
				vdm_down[j] += vdm[j*dmds+jj];
				vSdm_down[j] += vSdm[j*dmds+jj];
			}
			vdm_down[j] /= dmds;
			vSdm_down[j] /= dmds;
		}
	}

	/** plot */
	plt::Figure fig(8., 1.5);

	fig.set_background_color("black");
	fig.set_default_color("white");

	float adjustx = 0., adjusty = 0.02;
	/* profile */
	plt::Axes ax_pro(0.1+adjustx, 0.75+adjustx, 0.65+adjusty, 0.80+adjusty);
	ax_pro.plot(vt_down, vp_down);
	ax_pro.annotate("S/N = "+s_snr, 0.70, 0.8);
	ax_pro.annotate("w = "+s_width+" ms", 0.70, 0.65);
	ax_pro.axvline(vt[maxtid], 0., 1., {{"color", "red"}});
	ax_pro.autoscale(true, "x", true);
	ax_pro.set_ylabel("Flux");
	ax_pro.label(true, false, false, false);
	fig.push(ax_pro);

	/* f-t */
	plt::Axes ax_ft(0.1+adjustx, 0.75+adjustx, 0.35+adjusty, 0.65+adjusty);
	ax_ft.pcolor(vt_down, vf_down, mxft_down, "viridis");
	ax_ft.set_ylabel("Frequency (MHz)");
	ax_ft.label(true, false, false, false);
	fig.push(ax_ft);

	/* power-f */
	plt::Axes ax_powf(0.75+adjustx, 0.95+adjustx, 0.35+adjusty, 0.65+adjusty);
	ax_powf.plot(vpfcum_down, vfnorm_down);
	ax_powf.annotate("Integral Flux", 0.1, 1.2);
	ax_powf.axhline(fl, 0., 1., {{"color", "red"}});
	ax_powf.axhline(fh, 0., 1., {{"color", "red"}});
	ax_powf.autoscale(true, "y", true);
	ax_powf.label(false, false, false, true);
	fig.push(ax_powf);

	if (ndm > 1)
	{
		/* dm-t */
		float xpos = vt[maxtid];
		float ypos = vdm[maxdmid];
		float dmpos = dm_maxsnr;

		plt::Axes ax_dmt(0.1+adjustx, 0.75+adjustx, 0.05+adjusty, 0.35+adjusty);
		ax_dmt.pcolor(vt_down, vdm_down, mxdmt_down, "viridis");
		ax_dmt.circle(xpos, ypos, 5.);
		ax_dmt.set_xlabel("Time (s)");
		ax_dmt.set_ylabel("DM (cm\\u-3\\dpc)");
		fig.push(ax_dmt);

		/* dm-S */
		plt::Axes ax_dmS(0.75+adjustx, 0.95+adjustx, 0.05+adjusty, 0.35+adjusty);
		ax_dmS.plot(vSdm_down, vdm_down);
		float vSdm_down_min = *std::min_element(vSdm_down.begin(), vSdm_down.end());
		float vSdm_down_max = *std::max_element(vSdm_down.begin(), vSdm_down.end());
		float vdm_down_min = *std::min_element(vdm_down.begin(), vdm_down.end());
		ax_dmS.annotate("DM = "+s_dm+" cm\\u-3\\dpc", std::min(vSdm_down_min, threS)+0.75*(vSdm_down_max-std::min(vSdm_down_min, threS)), std::min(vdm_down_min, dmpos), {{"xycoords", "data"}, {"rotation", "270"}, {"refpos", "right"}});
		ax_dmS.axhline(dmpos, 0., 1., {{"color", "red"}});
		ax_dmS.axvline(threS, 0., 1., {{"color", "red"}});
		ax_dmS.set_xlabel("S/N");
		ax_dmS.label(false, false, true, false);
		ax_dmS.autoscale(true, "y", true);
		fig.push(ax_dmS);
	}

	/* metadata */
	plt::Axes ax_meta(0.1+adjustx, 0.95+adjustx, 0.81+adjusty, 0.99);
	ax_meta.label(false, false, false, false);
	ax_meta.frame(false, false, false, false);
	ax_meta.minorticks_off();
	ax_meta.majorticks_off();
	std::string fontsize = "1";
	ax_meta.annotate("Telescope = " + obsinfo["Telescope"], 0.02, 0.88, {{"fontsize", fontsize}});
	ax_meta.annotate("Beam = " + obsinfo["Beam"], 0.02, 0.74, {{"fontsize", fontsize}});
	ax_meta.annotate("RA = " + obsinfo["RA"], 0.02, 0.60, {{"fontsize", fontsize}});
	ax_meta.annotate("DEC = " + obsinfo["DEC"], 0.02, 0.46, {{"fontsize", fontsize}});
	ax_meta.annotate("DM (pc/cc) = " + s_dm, 0.02, 0.32, {{"fontsize", fontsize}});
	ax_meta.annotate("Width (ms) = " + s_width, 0.02, 0.18, {{"fontsize", fontsize}});
	ax_meta.annotate("S/N = " + s_snr, 0.02, 0.04, {{"fontsize", fontsize}});
	
	ax_meta.annotate("Source name = " + obsinfo["Source_name"], 0.42, 0.88, {{"fontsize", fontsize}});
	ax_meta.annotate("Date (MJD) = " + s_toa, 0.42, 0.74, {{"fontsize", fontsize}});
	ax_meta.annotate("GL (deg) = " + obsinfo["GL"], 0.42, 0.60, {{"fontsize", fontsize}});
	ax_meta.annotate("GB (deg) = " + obsinfo["GB"], 0.42, 0.46, {{"fontsize", fontsize}});
	ax_meta.annotate("MaxDM YMW16  (pc/cc) = " + s_ymw16_maxdm, 0.42, 0.32, {{"fontsize", fontsize}});
	ax_meta.annotate("Distance YMW16 (pc) = " + s_dist, 0.42, 0.18, {{"fontsize", fontsize}});
	
	ax_meta.annotate(rawdata, 0.985, 0.05, {{"xycoords","figure fraction"}, {"fontsize", "0.7"}, {"rotation", "270"}});
	ax_meta.annotate("Generated by TransientX and PlotX", 0.01, 0.01, {{"xycoords","figure fraction"}, {"fontsize", "0.7"}, {"rotation", "0"}});
	fig.push(ax_meta);

	fig.save(figname+"/PNG");

	std::ofstream outfile;
	outfile.open(rootname + "_replot.cands", ios_base::app); // append instead of overwrite
	
	outfile<<s_ibeam<<"\t";
	outfile<<s_k<<"\t";
	outfile<<setprecision(15)<<fixed<<mjd<<"\t";
	outfile<<setprecision(2)<<fixed<<dm_maxsnr<<"\t";
	outfile<<setprecision(2)<<fixed<<wn*tbin*1000<<"\t";
	outfile<<setprecision(1)<<fixed<<snr<<"\t";
	outfile<<setprecision(0)<<fixed<<fl<<"\t";
	outfile<<setprecision(0)<<fixed<<fh<<"\t";
	outfile<<figname<<"\t";
	outfile<<s_id<<"\t";
	outfile<<rawdata<<endl;
}

void Candidate::plot()
{
	/* prepare plot data */

	std::vector<float> vt(nbin, 0.);
	for (long int i=0; i<nbin; i++)
	{
		vt[i] = i*tbin;
	}

	std::vector<float> vf(nchan, 0.);
	for (long int j=0; j<nchan; j++)
	{
		vf[j] = frequencies[j];
	}

	std::vector<float> mxft(nchan*nbin, 0.);
	for (long int k=0; k<npol; k++)
	{
		if (k<2)
		{
			for (long int j=0; j<nchan; j++)
			{
				for (long int i=0; i<nbin; i++)
				{
					mxft[j*nbin+i] += data[k*nchan*nbin+j*nbin+i];
				}
			}
		}
	}

	std::ofstream outfile("tmp.dat", std::ios::binary);
	outfile.write((char *)mxft.data(), sizeof(float)*nchan*nbin);
	outfile.close();

	/* plot data */

	plt::Figure fig(8., 1.);

	fig.set_background_color("black");
	fig.set_default_color("white");

	float adjustx = 0., adjusty = 0.;

	/* f-t */
	plt::Axes ax_ft(0.1+adjustx, 0.95+adjustx, 0.1+adjusty, 0.95+adjusty);
	ax_ft.pcolor(vt, vf, mxft, "viridis");
	ax_ft.set_xlabel("Time (s)");
	ax_ft.set_ylabel("Frequency (MHz)");
	ax_ft.label(true, false, true, false);
	fig.push(ax_ft);

	fig.show();
}