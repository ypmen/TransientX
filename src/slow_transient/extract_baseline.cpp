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
#include "rescale.h"
#include "databuffer.h"

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
			("statsf", value<std::vector<std::string>>()->multitoken()->composing(), "stat list")
			("cont", "Input files are contiguous")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
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
			reader->sumif = false;
			reader->contiguous = contiguous;
			reader->verbose = verbose;
			reader->apply_scloffs = apply_scloffs;
			reader->apply_wts = apply_wts;
			reader->apply_zero_off = apply_zero_off;
			reader->check();
			reader->read_header();

			readers.push_back(std::move(reader));
		}
		else
		{
			auto reader = std::unique_ptr<FilterbankReader>(new FilterbankReader);
			reader->fnames = fnames;
			reader->sumif = false;
			reader->contiguous = contiguous;
			reader->verbose = verbose;
			reader->apply_scloffs = apply_scloffs;
			reader->apply_wts = apply_wts;
			reader->apply_zero_off = apply_zero_off;
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

	// read stats
	
	std::vector<HistogramAnalyzer> has;

	std::vector<float> chmean;
	std::vector<float> chstd;
	std::vector<float> chweight;

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
	else if (vm.count("statsf"))
	{
		std::vector<std::string> stats = vm["statsf"].as<std::vector<std::string>>();

		int nfile = stats.size();

		chweight.resize(nfile * nifs * nchans, 0.);
		chmean.resize(nfile * nifs * nchans, 0.);
		chstd.resize(nfile * nifs * nchans, 0.);

		for (auto stat=stats.begin(); stat!=stats.end(); ++stat)
		{
			std::ifstream fstat(*stat, std::ios::binary);

			fstat.read((char *)chweight.data(), chweight.size() * sizeof(float));
			fstat.read((char *)chmean.data(), chmean.size() * sizeof(float));
			fstat.read((char *)chstd.data(), chstd.size() * sizeof(float));
			fstat.close();
		}
	}

	std::vector<float> zerodm_x;
	std::vector<float> zerodm_y;
	if (nifs == 1)
	{
		zerodm_x.resize(ndump, 0.);
	}
	else
	{
		zerodm_x.resize(ndump, 0.);
		zerodm_y.resize(ndump, 0.);
	}

	DataBuffer<float> databuf(ndump, nifs * nchans);
	databuf.tsamp = tsamp;
	databuf.frequencies = readers[0]->frequencies;

	std::vector<Rescale> rescales;
	for (size_t k=0; k<readers.size(); k++)
	{
		Rescale rescale;
		rescale.prepare(databuf);

		std::copy(chmean.begin(), chmean.begin() + (k + 1) * nifs * nchans, rescale.chmean.begin());
		std::copy(chstd.begin(), chstd.begin() + (k + 1) * nifs * nchans, rescale.chstd.begin());
		std::copy(chweight.begin(), chweight.begin() + (k + 1) * nifs * nchans, rescale.chweight.begin());

		rescales.push_back(rescale);
	}

	Filterbank fil_x, fil_y;
	
	if (nifs == 1)
	{
		fil_x.filename = rootname + "_0.fil";
		std::strcpy(fil_x.source_name, source_name.c_str());
		fil_x.tstart = tstart;
		fil_x.tsamp = tsamp;
		fil_x.nchans = 1;
		fil_x.nbits = 32;
		fil_x.nifs = nifs;
		fil_x.fch1 = fch1;
		fil_x.foff = foff;
		fil_x.data_type = 2;
		fil_x.refdm = 0.;
		fil_x.refrm = 0.;

		if (!fil_x.write_header())
		BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";
	}
	else
	{
		fil_x.filename = rootname + "_0.fil";
		std::strcpy(fil_x.source_name, source_name.c_str());
		fil_x.tstart = tstart;
		fil_x.tsamp = tsamp;
		fil_x.nchans = 1;
		fil_x.nbits = 32;
		fil_x.nifs = 1;
		fil_x.fch1 = fch1;
		fil_x.foff = foff;
		fil_x.data_type = 2;
		fil_x.refdm = 0.;
		fil_x.refrm = 0.;

		fil_y.filename = rootname + "_1.fil";
		std::strcpy(fil_y.source_name, source_name.c_str());
		fil_y.tstart = tstart;
		fil_y.tsamp = tsamp;
		fil_y.nchans = 1;
		fil_y.nbits = 32;
		fil_y.nifs = 1;
		fil_y.fch1 = fch1;
		fil_y.foff = foff;
		fil_y.data_type = 2;
		fil_y.refdm = 0.;
		fil_y.refrm = 0.;

		if (!fil_x.write_header())
		BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";

		if (!fil_y.write_header())
		BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";
	}
	
	while (true)
	{
		if (nifs == 1)
		{
			std::fill(zerodm_x.begin(), zerodm_x.end(), 0.);
		}
		else
		{
			std::fill(zerodm_x.begin(), zerodm_x.end(), 0.);
			std::fill(zerodm_y.begin(), zerodm_y.end(), 0.);
		}

		int id = 0;
		for (auto reader=readers.begin(); reader!=readers.end(); ++reader)
		{
			if ((*reader)->is_end) goto end;

			if ((*reader)->read_data(databuf, ndump) != ndump) goto end;

			DataBuffer<float> *data = NULL;

			if (vm.count("stats"))
			{
				has[id].filter(databuf);

				data = databuf.get();
			}
			else if (vm.count("statsf"))
			{
				data = rescales[id].filter(databuf);
			}
			else
			{
				data = databuf.get();
			}

			if (nifs == 1)
			{
				for (size_t i=0; i<ndump; i++)
				{
					for (size_t j=0; j<nchans; j++)
					{
						zerodm_x[i] += data->buffer[i * nchans + j];
					}
				}
			}
			else
			{
				for (size_t i=0; i<ndump; i++)
				{
					for (size_t j=0; j<nchans; j++)
					{
						zerodm_x[i] += data->buffer[i * nifs* nchans + 0 * nchans + j];
						zerodm_y[i] += data->buffer[i * nifs* nchans + 1 * nchans + j];
					}
				}
			}

			id++;
		}

		if (nifs == 1)
		{
			fwrite(zerodm_x.data(), 1, sizeof(float) * ndump, fil_x.fptr);
		}
		else
		{
			fwrite(zerodm_x.data(), 1, sizeof(float) * ndump, fil_x.fptr);
			fwrite(zerodm_y.data(), 1, sizeof(float) * ndump, fil_y.fptr);
		}
	}
	end:

	if (nifs == 1)
	{
		fil_x.close();
	}
	else
	{
		fil_x.close();
		fil_y.close();
	}
}