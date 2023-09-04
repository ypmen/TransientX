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
			("stat", value<std::string>()->multitoken()->composing(), "histo list")
			("statf", value<std::string>()->multitoken()->composing(), "stat list")
			("td", value<int>()->default_value(1), "Time downsample")
			("fd", value<int>()->default_value(1), "Frequency downsample")
			("DMs", value<std::vector<double>>()->multitoken()->composing(), "Dedispersion measure")
			("RMs", value<std::vector<double>>()->multitoken()->composing(), "Rotation measure")
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

	int td = vm["td"].as<int>();
	int fd = vm["fd"].as<int>();

	std::string rootname = vm["output"].as<std::string>();

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

	int ndump = 4096;

	assert(ndump % td == 0);
	assert(nchans % fd == 0);

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

	Rescale rescale;
	rescale.prepare(databuf);
	std::copy(chmean.begin(), chmean.end(), rescale.chmean.begin());
	std::copy(chstd.begin(), chstd.end(), rescale.chstd.begin());
	std::copy(chweight.begin(), chweight.end(), rescale.chweight.begin());

	std::vector<Defaraday> defaradays;
	std::vector<Dedispersion> dedispersions;

	for (size_t k=0; k<DMs.size(); k++)
	{
		Defaraday defaraday;
		defaraday.rm = RMs[k];
		defaraday.prepare(rescale);

		Dedispersion dedispersion;
		dedispersion.dm = DMs[k];
		dedispersion.prepare(defaraday);

		defaradays.push_back(defaraday);
		dedispersions.push_back(dedispersion);
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
		fils[k].tsamp = tsamp * td;
		fils[k].nchans = nchans / fd;
		fils[k].nbits = 32;
		fils[k].nifs = nifs;
		fils[k].fch1 = fch1;
		fils[k].foff = foff;
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
			DataBuffer<float> *data = defaradays[k].run(databuf);

			data = dedispersions[k].run(*data);

			std::vector<float> temp((ndump / td) * nifs * (nchans / fd), 0.);

			size_t ndump_td = ndump / td;
			size_t nchans_fd = nchans / fd;

			if (data->counter > dedispersions[k].offset)
			{
				for (size_t i=0; i<ndump_td; i++)
				{
					for (size_t n=0; n<td; n++)
					{
						for (size_t l=0; l<nifs; l++)
						{
							for (size_t j=0; j<nchans_fd; j++)
							{
								for (size_t m=0; m<fd; m++)
								{
									temp[i * nifs * nchans_fd + l * nchans_fd + j] += data->buffer[(i * td + n) * nifs * nchans + l * nchans + (j * fd + m)];
								}
							}
						}
					}
				}

				fwrite(temp.data(), 1, sizeof(float) * temp.size(), fils[k].fptr);
			}
		}
	}

	for (size_t k=0; k<DMs.size(); k++) fils[k].close();
}