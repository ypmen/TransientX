/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2023-07-10 11:29:38
 * @modify date 2023-07-10 11:29:38
 * @desc [description]
 */

#include <iostream>
#include <fstream>
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
			("stats", value<std::vector<std::string>>()->multitoken()->composing(), "histo list")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
			("output,o", value<std::string>()->default_value("baseline"), "Output rootname")
			("list,l", value<std::string>(), "Input list file");

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
	if (vm.count("list") == 0)
	{
		std::cerr<<"Error: no input file"<<std::endl;
		return -1;
	}

	bool contiguous = vm.count("cont");

	num_threads = vm["threads"].as<unsigned int>();

	std::string rootname = vm["output"].as<std::string>();

	// read stats
	
	std::vector<HistogramAnalyzer> has;

	if (vm.count("stats"))
	{
		std::vector<std::string> stats = vm["stats"].as<std::vector<std::string>>();

		for (auto stat=stats.begin(); stat!=stats.end(); ++stat)
		{
			HistogramAnalyzer ha;
			ha.read(*stat);
			ha.get_statitstics();

			has.push_back(ha);
		}
	}

	// read data files
	std::string input_list = vm["list"].as<std::string>();

	std::vector<std::unique_ptr<PSRDataReader>> readers;

	std::ifstream input_list_f(input_list);
	std::string line;
	while (getline(input_list_f, line))
	{
		boost::trim(line);
		std::vector<string> fnames;
		boost::split(fnames, line, boost::is_any_of("\t "), boost::token_compress_on);

		if (vm.count("psrfits"))
		{
			auto reader = std::unique_ptr<PsrfitsReader>(new PsrfitsReader);
			reader->fnames = fnames;
			reader->sumif = true;
			reader->contiguous = contiguous;
			reader->verbose = verbose;
			reader->apply_scloffs = false;
			reader->apply_wts = false;
			reader->apply_zero_off = false;
			reader->check();
			reader->read_header();

			readers.push_back(std::move(reader));
		}
		else
		{
			auto reader = std::unique_ptr<FilterbankReader>(new FilterbankReader);
			reader->fnames = fnames;
			reader->sumif = true;
			reader->contiguous = contiguous;
			reader->verbose = verbose;
			reader->apply_scloffs = false;
			reader->apply_wts = false;
			reader->apply_zero_off = false;
			reader->check();
			reader->read_header();

			readers.push_back(std::move(reader));
		}
	}

	// align start sample
	long double tstart = 0.;
	for (auto reader=readers.begin(); reader!=readers.end(); ++reader)
	{
		assert((*reader)->tsamp == readers[0]->tsamp);
		assert((*reader)->nchans == readers[0]->nchans);

		long double t = (*reader)->start_mjd.to_day();
		if (t > tstart)
			tstart = t;
	}

	for (auto reader=readers.begin(); reader!=readers.end(); ++reader)
	{
		long double t = (*reader)->start_mjd.to_day();
		(*reader)->skip_start = std::round((tstart - t) * 86400. / (*reader)->tsamp);
		(*reader)->skip_head();
	}

	size_t nchans = readers[0]->nchans;
	double tsamp = readers[0]->tsamp;
	int nifs = readers[0]->nifs;
	double fch1 = 0.;
	double foff = 0.;
	readers[0]->get_fch1_foff(fch1, foff);
	std::string source_name = readers[0]->source_name;

	int ndump = 1024;

	std::vector<float> zerodm(ndump, 0.);
	DataBuffer<float> databuf(ndump, nchans);

	Filterbank fil;
	fil.filename = rootname + ".fil";
	std::strcpy(fil.source_name, source_name.c_str());
	fil.tstart = tstart;
	fil.tsamp = tsamp;
	fil.nchans = 1;
	fil.nbits = 32;
	fil.nifs = nifs;
	fil.fch1 = fch1;
	fil.foff = foff;
	fil.data_type = 2;
	fil.refdm = 0.;

	if (!fil.write_header())
		BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";

	while (true)
	{
		std::fill(zerodm.begin(), zerodm.end(), 0.);

		int id = 0;
		for (auto reader=readers.begin(); reader!=readers.end(); ++reader)
		{
			if ((*reader)->is_end) goto end;

			if ((*reader)->read_data(databuf, ndump) != ndump) goto end;

			if (vm.count("stats"))
			{
				has[id].filter(databuf);
			}

			for (size_t i=0; i<ndump; i++)
			{
				for (size_t j=0; j<nchans; j++)
				{
					zerodm[i] += databuf.buffer[i * nchans + j];
				}
			}

			id++;
		}

		fwrite(zerodm.data(), 1, sizeof(float) * ndump, fil.fptr);
	}
	end:

	fil.close();
}