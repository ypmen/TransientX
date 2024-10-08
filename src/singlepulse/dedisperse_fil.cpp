/*
 * dedisperse.cpp
 *
 *  Created on: Apr 24, 2020
 *      Author: ypmen
 */

#define EIGHTBIT 1
#define NSBLK 65536

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <iomanip>
#include <boost/program_options.hpp>
#include <random>

#include "dedisperse.h"
#include "psrfits.h"
#include "integration.h"
#include "mjd.h"
#include "databuffer.h"
#include "singlepulse.h"
#include "patch.h"
#include "preprocesslite.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "logging.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
unsigned int dbscan_radius;
unsigned int dbscan_k;

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
			("td", value<int>()->default_value(1), "Time downsample for preprocessing")
			("fd", value<int>()->default_value(1), "Frequency downsample for preprocessing")
			("zapthre", value<float>()->default_value(4), "Threshold in IQR for zapping channels")
			("dms", value<double>()->default_value(0), "DM start")
			("ddm", value<double>()->default_value(1), "DM step")
			("ndm", value<int>()->default_value(200), "Number of DM")
			("overlap", value<double>()->default_value(0), "Overlap ratio")
			("ddplan", value<string>(), "Input ddplan file")
			("thre", value<float>()->default_value(7), "S/N threshold")
			("minw", value<float>()->default_value(1e-4), "Minimum pulse width (s)")
			("maxw", value<float>()->default_value(0.02), "Maximum pulse width (s)")
			("snrloss", value<float>()->default_value(1e-1), "S/N loss")
			("seglen,l", value<float>()->default_value(1), "Time length per segment (s)")
			("ra", value<string>()->default_value("00:00:00"), "RA (hh:mm:ss.s)")
			("dec", value<string>()->default_value("00:00:00"), "DEC (dd:mm:ss.s)")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("telescope", value<string>()->default_value("Fake"), "Telescope name")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("radius,r", value<double>()->default_value(1), "DBSCAN radius (ms)")
			("neighbors,k", value<unsigned int>()->default_value(2), "DBSCAN k")
			("maxncand", value<int>()->default_value(100), "Maximum number of candidates in one data block")
			("minpts", value<int>()->default_value(5), "Minimum points in one cluster")
			("baseline", value<vector<float>>()->multitoken()->default_value(vector<float>{0.0, 0.1}, "0.0, 0.1"), "The scale of baseline remove (s)")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threPatch", value<float>()->default_value(10), "IQR threshold of patch for bad data")
			("widthPatch", value<float>()->default_value(0.4), "Width threshold (s) of patch for bad data")
			("fillPatch", value<std::string>()->default_value("none"), "Fill the bad data by [none, mean, rand] in patch")
			("fill", value<string>()->default_value("mean"), "Fill the zapped samples by [mean, rand]")
			("source_name,s", value<string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("drop", "Drop candidates with maximum search width")
			("iqr", "Calculate variance and mean based on IQR")
			("saveimage", "Save images to fits")
			("repeater", "Using 2D matched filter (under development)")
			("mean", value<float>()->default_value(0), "Mean of dedispersed time series")
			("std", value<float>()->default_value(3), "Standard deviation of dedispersed time series")
			("nbits", value<int>()->default_value(8), "Data type of dedispersed time series")
			("savetim", "Output dedispersed data (sigproc format by default)")
			("format", value<string>()->default_value("pulsarx"), "Output format of dedispersed data [pulsarx(default),sigproc,presto]")
			("config,c", value<std::string>(), "Config json file (This will cover all relevant configurations)")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
			("input,f", value<vector<string>>()->multitoken()->composing(), "Input files");

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

	nlohmann::json config;
	nlohmann::json config_patch;
	nlohmann::json config_prep;

	if (vm.count("config"))
	{
		std::ifstream config_f(vm["config"].as<std::string>());
		config = nlohmann::json::parse(config_f);

		config_patch = config["patch"];
		config_prep = config["preprocesslite"];
	}

	bool contiguous = vm.count("cont");

	num_threads = vm["threads"].as<unsigned int>();

	vector<double> jump = vm["jump"].as<vector<double>>();

	vector<string> fnames = vm["input"].as<vector<string>>();

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

	long double tstart = reader->start_mjd.to_day();

	string source_name = reader->source_name;
	string s_telescope = reader->telescope;
	int ibeam = reader->beam.empty() ? 0 : std::stoi(reader->beam);
	double src_raj = 0., src_dej;
	if (!reader->ra.empty())
	{
		string ra = reader->ra;
		ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
		src_raj = stod(ra);
	}

	if (!reader->dec.empty())
	{
		string dec = reader->dec;
		dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
		src_dej = stod(dec);
	}

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	vector<SinglePulse> search1;
	if (vm.count("config"))
		parse_json(vm, config, search1);
	else
		parse(vm, search1);

	vector<int> tds;
	if (vm.count("config"))
	{
		for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
		{
			tds.push_back((*sp).td*(long int)(config_prep["td"]));
		}
	}
	else
	{
		for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
		{
			tds.push_back((*sp).td*vm["td"].as<int>());
		}
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm*td_lcm;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader->frequencies.data(), sizeof(double)*nchans);

	Patch patch;
	if (vm.count("config"))
	{
		patch.read_config(config_patch);
	}
	else
	{
		patch.filltype = vm["fillPatch"].as<string>();
		patch.width = vm["widthPatch"].as<float>();
		patch.threshold = vm["threPatch"].as<float>();
	}
	patch.prepare(databuf);
	patch.close();

	PreprocessLite prep;
	if (vm.count("config"))
	{
		prep.read_config(config_prep);
	}
	else
	{
		prep.td = vm["td"].as<int>();
		prep.fd = vm["fd"].as<int>();
		prep.thresig = vm["zapthre"].as<float>();
		prep.filltype = vm["fill"].as<string>();
	}
	prep.prepare(databuf);

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	reader->skip_start = nstart;
	reader->skip_end = ntotal - nend;
	reader->skip_head();

	stringstream ss_ibeam;
	if (vm.count("incoherent"))
		ss_ibeam << "ifbf" << setw(5) << setfill('0') << ibeam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << ibeam;
	string s_ibeam = ss_ibeam.str();

	long int nsearch = search1.size();
	for (long int k=0; k<nsearch; k++)
	{
		search1[k].tstart = tstart+nstart*tsamp/86400.;
		search1[k].source_name = source_name;
		search1[k].telescope = s_telescope;
		search1[k].ibeam = ibeam;
		search1[k].src_raj = src_raj;
		search1[k].src_dej = src_dej;

		reader->get_filterbank_template(search1[k].fildedisp);
		search1[k].fildedisp.tstart = tstart+nstart*tsamp/86400.;
		search1[k].fildedisp.ibeam = ibeam;
		search1[k].fildedisp.fch1 = prep.frequencies.front();
		search1[k].fildedisp.foff = prep.frequencies.back()-databuf.frequencies.front();
		search1[k].fildedisp.nchans = prep.frequencies.size();
		search1[k].prepare(prep);
	}

	if (verbose == 1)
	{
		cout<<"Maximum width = "<<setprecision(1)<<fixed<<search1[0].maxw*1000<<" (ms)"<<endl;
	}

	while (!reader->is_end)
	{
		long int idxn = reader->get_ifile();
		long int n = reader->get_ifile_ordered();

		if (reader->read_data(databuf, ndump) != ndump) break;

		databuf.counter += ndump;

		patch.filter(databuf);
		prep.run(databuf);
		for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
		{
			prep.isbusy = true;
			(*sp).fileid = idxn+1;
			(*sp).fname = fnames[n];
			(*sp).verbose = verbose;
			(*sp).run(prep);
		}
	}

	std::random_device rd{};
    std::mt19937 gen{rd()};
	std::normal_distribution<> distribution{0, 1};
	std::generate(prep.buffer.begin(), prep.buffer.end(), [&distribution, &gen](){return distribution(gen);});

	for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
	{
		int nleft = sp->dedisp.offset/sp->dedisp.ndump+1;
		for (long int k=0; k<nleft; k++)
		{
			prep.counter += prep.nsamples;
			prep.isbusy = true;
			sp->run(prep);
		}
	}

	if (vm["format"].as<string>() == "presto")
	{
		for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
		{
			(*sp).dedisp.makeinf((*sp).fildedisp);
		}
	}
	else if (vm["format"].as<string>() != "sigproc")
	{
		for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
		{
			(*sp).dedisp.modifynblock();
		}
	}

	return 0;
}
