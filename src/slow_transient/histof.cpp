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
#include "stat.h"

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
			("zapthre", value<float>()->default_value(3), "Threshold in IQR for zapping channels")
			("nosumif", "Keep polarization")
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

	bool nosumif = vm.count("nosumif");

	std::vector<double> jump = vm["jump"].as<std::vector<double>>();

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = !nosumif;
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

	long int nstart = jump[0]/tsamp;
	long int nend = jump[1]/tsamp;

	reader->skip_start = nstart;
	reader->skip_end = nend;
	reader->skip_head();

	tstart += (nstart * tsamp) / 86400.;

	int ndump = 1024;

	DataBuffer<float> databuf;
	if (nosumif)
		databuf.resize(ndump, nifs * nchans);
	else
		databuf.resize(ndump, nchans);

	Stat stat;
	stat.zap_threshold = vm["zapthre"].as<float>();
	stat.prepare(databuf);

	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		stat.run(databuf);
	}

	stat.get_stat();

	std::ofstream outfile;
	outfile.open(rootname + ".dat", std::ios::binary);
	outfile.write((char *)stat.chweight.data(), stat.chweight.size() * sizeof(float));
	outfile.write((char *)stat.chmean.data(), stat.chmean.size() * sizeof(float));
	outfile.write((char *)stat.chstd.data(), stat.chstd.size() * sizeof(float));
	outfile.close();

	delete reader;

	return 0;
}