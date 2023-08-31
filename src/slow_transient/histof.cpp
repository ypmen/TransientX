/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2023-06-21 20:19:04
 * @modify date 2023-06-21 20:19:04
 * @desc [description]
 */

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
using namespace boost::program_options;

#include "logging.h"
#include "psrdatareader.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("zapthre", value<float>()->default_value(3), "Threshold in IQR for zapping channels")
			("cont", "Input files are contiguous")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("psrfits", "Input psrfits format data")
			("output,o", value<std::string>()->default_value("stat"), "Output rootname")
			("input,f", value<std::vector<std::string>>()->multitoken()->composing(), "Input files");

	positional_options_description pos_desc;
	pos_desc.add("input", -1);
	command_line_parser parser{argc, argv};
	parser.options(desc).style(command_line_style::default_style | command_line_style::allow_short);
	parser.options(desc).positional(pos_desc);
	parsed_options parsed_options = parser.run();

	variables_map vm;
	store(parsed_options, vm);
	notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << '\n';
		return 0;
	}
	if (vm.count("verbose"))
	{
		verbose = 1;
	}
	if (vm.count("input") == 0)
	{
		std::cerr<<"Error: no input file"<<std::endl;
		return -1;
	}

	bool contiguous = vm.count("cont");

	bool apply_wts = false;
	bool apply_scloffs = false;
	bool apply_zero_off = false;

	if (vm.count("wts"))
		apply_wts = true;
	if (vm.count("scloffs"))
		apply_scloffs = true;
	if (vm.count("zero_off"))
		apply_zero_off = true;

	num_threads = vm["threads"].as<unsigned int>();

	std::string rootname = vm["output"].as<std::string>();

	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	long double tstart = reader->start_mjd.to_day();
	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	int ndump = 1024;

	DataBuffer<float> databuf(ndump, nchans);

	std::vector<float> chmean(nchans, 0.);
	std::vector<float> chstd(nchans, 0.);
	std::vector<float> chskewness(nchans, 0.);
	std::vector<float> chkurtosis(nchans, 0.);
	std::vector<float> chweight(nchans, 0.);

	std::vector<double> chmean1(nchans, 0.);
	std::vector<double> chmean2(nchans, 0.);
	std::vector<double> chmean3(nchans, 0.);
	std::vector<double> chmean4(nchans, 0.);
	std::vector<double> chcorr(nchans, 0.);
	std::vector<double> last_data(nchans, 0.);

	size_t count = 0;
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		for (size_t i=0; i<databuf.nsamples; i++)
		{
			for (size_t j=0; j<databuf.nchans; j++)
			{
				double tmp1 = databuf.buffer[i * databuf.nchans + j];
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

		count += ndump;
	}

	for (long int j=0; j<databuf.nchans; j++)
	{
		chmean1[j] /= count;
		chmean2[j] /= count;
		chmean3[j] /= count;
		chmean4[j] /= count;

		chcorr[j] /= count - 1;

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
	std::vector<float> kurtosis_sort(chkurtosis.begin(), chkurtosis.end());
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::less<float>());
	float kurtosis_q1 = kurtosis_sort[kurtosis_sort.size()/4];
	std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::greater<float>());
	float kurtosis_q3 =kurtosis_sort[kurtosis_sort.size()/4];
	float kurtosis_R = kurtosis_q3-kurtosis_q1;

	std::vector<float> skewness_sort(chskewness.begin(), chskewness.end());
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

	float thresig = vm["zapthre"].as<float>();

	for (long int j=0; j<nchans; j++)
	{
		if (chkurtosis[j]>=kurtosis_q1-thresig*kurtosis_R && \
			chkurtosis[j]<=kurtosis_q3+thresig*kurtosis_R && \
			chskewness[j]>=skewness_q1-thresig*skewness_R && \
			chskewness[j]<=skewness_q3+thresig*skewness_R && \
			chcorr[j]>=corr_q1-thresig*corr_R && \
			chcorr[j]<=corr_q3+thresig*corr_R)
		{
			chweight[j] = 1.;
		}
	}

	std::ofstream outfile;
	outfile.open(rootname + ".dat", std::ios::binary);
	outfile.write((char *)chweight.data(), chweight.size() * sizeof(float));
	outfile.write((char *)chmean.data(), chmean.size() * sizeof(float));
	outfile.write((char *)chstd.data(), chstd.size() * sizeof(float));
	outfile.close();

	delete reader;

	return 0;
}