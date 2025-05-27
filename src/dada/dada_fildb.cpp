/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-05 09:00:52
 * @modify date 2024-04-05 09:00:52
 * @desc "read filterbank data to psrdada buffer"
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "dada.h"
#include "databuffer.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "logging.h"

using namespace boost::program_options;

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
			("jump,j", value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("ra", value<std::string>()->default_value("00:00:00"), "RA (hh:mm:ss.s)")
			("dec", value<std::string>()->default_value("00:00:00"), "DEC (dd:mm:ss.s)")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("source_name,s", value<std::string>()->default_value("J0000-00"), "Source name")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
			("key_output,k", value<std::string>()->default_value("1111"), "Output dada key")
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
		cerr<<"Error: no input file"<<endl;
		return -1;
	}

	bool contiguous = vm.count("cont");

	num_threads = vm["threads"].as<unsigned int>();

	std::vector<double> jump = vm["jump"].as<std::vector<double>>();

	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	bool apply_wts = false;
	bool apply_scloffs = false;
	bool apply_zero_off = false;

	if (vm.count("wts"))
		apply_wts = true;
	if (vm.count("scloffs"))
		apply_scloffs = true;
	if (vm.count("zero_off"))
		apply_zero_off = true;

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	if (!vm["ibeam"].defaulted())
		reader->beam = std::to_string(vm["ibeam"].as<int>());
	if (!vm["telescope"].defaulted())
		reader->telescope = vm["telescope"].as<std::string>();
	if (!vm["ra"].defaulted())
		reader->ra = vm["ra"].as<std::string>();
	if (!vm["dec"].defaulted())
		reader->dec = vm["dec"].as<std::string>();
	if (!vm["source_name"].defaulted())
		reader->source_name = vm["source_name"].as<std::string>();

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	reader->skip_start = nstart;
	reader->skip_end = ntotal - nend;
	reader->skip_head();

	// write output header
	stringstream ss_tstart;
	ss_tstart << setprecision(13) << fixed << reader->start_mjd.to_day();
	string s_tstart = ss_tstart.str();

	nlohmann::json output_header = {
		{"telescope", reader->telescope},
		{"source_name",reader->source_name},
		{"ra", reader->ra},
		{"dec", reader->dec},
		{"beam", reader->beam},
		{"tstart", s_tstart},
		{"nifs", reader->nifs},
		{"nbits", 8},
		{"nchans", reader->nchans},
		{"tsamp", reader->tsamp},
		{"fch1", reader->frequencies.front()},
		{"foff", reader->frequencies[1]-reader->frequencies[0]}
	};

	std::string key_output = vm["key_output"].as<std::string>();
	PSRDADA::Writer writer(key_output);
	writer.prepare(output_header);

	// create databuffer
	long int ndump = writer.get_bufsz() / nchans;

	DataBuffer<unsigned char> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(databuf.frequencies.data(), reader->frequencies.data(), sizeof(double)*nchans);

	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		databuf.counter += ndump;

		writer.run((char *)(databuf.buffer.data()), databuf.buffer.size());
	}
}