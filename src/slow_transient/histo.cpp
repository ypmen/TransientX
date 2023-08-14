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
			("cont", "Input files are contiguous")
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
	reader->apply_scloffs = false;
	reader->apply_wts = false;
	reader->apply_zero_off = false;
	reader->check();
	reader->read_header();

	long double tstart = reader->start_mjd.to_day();
	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	int ndump = 1024;

	DataBuffer<unsigned char> databuf(ndump, nchans);

	std::vector<size_t> hist(256 * nchans);
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		for (size_t i=0; i<databuf.nsamples; i++)
		{
			for (size_t j=0; j<databuf.nchans; j++)
			{
				unsigned char d = databuf.buffer[i * databuf.nchans + j];
				hist[d * nchans + j] += 1;
			}
		}
	}

	std::ofstream outfile;
	outfile.open(rootname + ".dat", std::ios::binary);
	outfile.write((char *)hist.data(), hist.size() * sizeof(size_t));
	outfile.close();

	delete reader;

	return 0;
}