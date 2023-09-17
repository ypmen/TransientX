/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2023-09-03 16:48:35
 * @modify date 2023-09-03 16:48:35
 * @desc [description]
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
using namespace boost::program_options;
#include <boost/algorithm/string.hpp>
#include <memory>

#include "logging.h"
#include "psrdatareader.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "histogram_analyzer.h"
#include "rescale.h"
#include "defaraday.h"
#include "dedispersion.h"
#include "downsample.h"

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
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("range", value<vector<double>>()->multitoken(), "start and end time (s)")
			("stat", value<std::string>()->multitoken()->composing(), "histo list")
			("statf", value<std::string>()->multitoken()->composing(), "stat list")
			("tds", value<std::vector<int>>()->multitoken()->composing(), "Time downsamples")
			("fds", value<std::vector<int>>()->multitoken()->composing(), "Frequency downsamples")
			("DMs", value<std::vector<double>>()->multitoken()->composing(), "Dedispersion measures")
			("RMs", value<std::vector<double>>()->multitoken()->composing(), "Rotation measures")
			("cont", "Input files are contiguous")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("psrfits", "Input psrfits format data")
			("output,o", value<std::string>()->default_value("source"), "Output rootname")
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

	std::vector<double> DMs;
	if (vm.count("DMs"))
	{
		DMs = vm["DMs"].as<std::vector<double>>();
	}
	else
	{
		DMs.push_back(0.);
	}

	std::vector<double> RMs;
	if (vm.count("RMs"))
	{
		RMs = vm["RMs"].as<std::vector<double>>();
	}
	else
	{
		RMs.push_back(0.);
	}

	std::vector<int> tds;
	if (vm.count("tds"))
	{
		tds = vm["tds"].as<std::vector<int>>();
	}
	else
	{
		tds.push_back(1);
	}

	std::vector<int> fds;
	if (vm.count("fds"))
	{
		fds = vm["fds"].as<std::vector<int>>();
	}
	else
	{
		fds.push_back(1);
	}

	assert(DMs.size() == RMs.size());

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

	std::vector<double> jump = vm["jump"].as<std::vector<double>>();

	// read data files
	vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = false;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	long double tstart = reader->start_mjd.to_day();
	size_t nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	double fch1 = 0.;
	double foff = 0.;
	reader->get_fch1_foff(fch1, foff);
	std::string source_name = reader->source_name;

	long int nstart = jump[0]/tsamp;
	long int nend = jump[1]/tsamp;
	
	if (vm.count("range"))
	{
		std::vector<double> range = vm["range"].as<std::vector<double>>();
		nstart = range[0] / tsamp;
		nend = reader->nsamples - (int)(range[1] / tsamp);
		std::cout<<nstart<<" "<<nend<<std::endl;
	}

	reader->skip_start = nstart;
	reader->skip_end = nend;
	reader->skip_head();

	tstart += (nstart * tsamp) / 86400.;

	int ndump = 4096;

	// read stats
	
	HistogramAnalyzer ha;

	std::vector<float> chmean;
	std::vector<float> chstd;
	std::vector<float> chweight;

	if (vm.count("stat"))
	{
		std::string stat = vm["stat"].as<std::string>();

		ha.read(stat);
		ha.get_statitstics();
	}
	else if (vm.count("statf"))
	{
		std::string stat = vm["statf"].as<std::string>();

		chweight.resize(nifs * nchans, 0.);
		chmean.resize(nifs * nchans, 0.);
		chstd.resize(nifs * nchans, 0.);

		std::ifstream fstat(stat, std::ios::binary);

		fstat.read((char *)chweight.data(), chweight.size() * sizeof(float));
		fstat.read((char *)chmean.data(), chmean.size() * sizeof(float));
		fstat.read((char *)chstd.data(), chstd.size() * sizeof(float));

		fstat.close();
	}

	// initalize kernels

	DataBuffer<float> databuf;
	databuf.resize(ndump, nifs * nchans);
	databuf.tsamp = tsamp;
	databuf.frequencies = reader->frequencies;
	databuf.means.resize(nifs * nchans, 0.);
	databuf.vars.resize(nifs * nchans, 0.);
	databuf.weights.resize(nifs * nchans, 0.);

	Rescale rescale;
	rescale.prepare(databuf);
	std::copy(chmean.begin(), chmean.end(), rescale.chmean.begin());
	std::copy(chstd.begin(), chstd.end(), rescale.chstd.begin());
	std::copy(chweight.begin(), chweight.end(), rescale.chweight.begin());

	std::vector<Downsample> tdownsamples;
	std::vector<Defaraday> defaradays;
	std::vector<Dedispersion> dedispersions;
	std::vector<Downsample> fdownsamples;

	for (size_t k=0; k<DMs.size(); k++)
	{
		assert(ndump % tds[k] == 0);
		assert(nchans % fds[k] == 0);

		Downsample tdownsample;
		tdownsample.td = tds[k];
		tdownsample.fd = 1;
		tdownsample.prepare(rescale);

		Defaraday defaraday;
		defaraday.rm = RMs[k];
		defaraday.prepare(tdownsample);

		Dedispersion dedispersion;
		dedispersion.dm = DMs[k];
		dedispersion.prepare(defaraday);

		Downsample fdownsample;
		fdownsample.td = 1;
		fdownsample.fd = fds[k];
		fdownsample.prepare(dedispersion);

		tdownsamples.push_back(tdownsample);
		defaradays.push_back(defaraday);
		dedispersions.push_back(dedispersion);
		fdownsamples.push_back(fdownsample);	
	}

	std::vector<Filterbank> fils;
	fils.resize(DMs.size());
	for (size_t k=0; k<DMs.size(); k++)
	{
		std::stringstream ss_dm;
		ss_dm << "DM" << setprecision(2) << fixed << setfill('0') << DMs[k];
		std::string s_dm = ss_dm.str();

		std::stringstream ss_rm;
		ss_rm << "RM" << setprecision(2) << fixed << setfill('0') << RMs[k];
		std::string s_rm = ss_rm.str();

		fils[k].filename = rootname + "_" + s_dm + "_" + s_rm  + ".fil";
		std::strcpy(fils[k].source_name, source_name.c_str());
		fils[k].tstart = tstart;
		fils[k].tsamp = tsamp * tds[k];
		fils[k].nchans = nchans / fds[k];
		fils[k].nbits = 32;
		fils[k].nifs = nifs;
		fils[k].fch1 = fch1 + (fds[k] - 1) * foff * 0.5;
		fils[k].foff = foff * fds[k];
		fils[k].data_type = 1;
		fils[k].refdm = DMs[k];
		fils[k].refrm = RMs[k];

		if (!fils[k].write_header())
			BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";
	}

	while (!reader->is_end)
	{	
		if (reader->read_data(databuf, ndump) != ndump) break;
		databuf.counter += ndump;

		if (vm.count("stat"))
		{
			ha.filter(databuf);
		}
		else if (vm.count("statf"))
		{
			rescale.filter(databuf);
		}

		for (size_t k=0; k<DMs.size(); k++)
		{
			DataBuffer<float> *data = tdownsamples[k].run(databuf);

			data = defaradays[k].run(*data);

			data = dedispersions[k].run(*data);

			data = fdownsamples[k].run(*data);

			if (data->counter > dedispersions[k].offset)
			{
				fwrite(data->buffer.data(), 1, sizeof(float) * data->buffer.size(), fils[k].fptr);
			}
		}
	}

	for (size_t k=0; k<DMs.size(); k++) fils[k].close();
}