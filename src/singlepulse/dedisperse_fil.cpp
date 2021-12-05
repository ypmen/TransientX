/*
 * dedisperse.cpp
 *
 *  Created on: May 10, 2020
 *      Author: ypmen
 */

#define FAST 1

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <iomanip>
#include <boost/program_options.hpp>

#include "dedisperse.h"
#include "filterbank.h"
#include "mjd.h"
#include "utils.h"
#include "databuffer.h"
#include "singlepulse.h"
#include "preprocesslite.h"
#include "logging.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
unsigned int dbscan_radius;
unsigned int dbscan_k;

bool repeater;
bool dumptim=false;

#define NSBLK 65536

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
			("ra", value<double>()->default_value(0), "RA (hhmmss.s)")
			("dec", value<double>()->default_value(0), "DEC (ddmmss.s)")
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
			("fill", value<string>()->default_value("mean"), "Fill the zapped samples by [mean, rand]")
			("source_name,s", value<string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("drop", "Drop candidates with maximum search width")
			("iqr", "Calculate variance and mean based on IQR")
			("repeater", "Using 2D matched filter (under development)")
			("mean", value<float>()->default_value(0), "Mean of dedispersed time series")
			("std", value<float>()->default_value(3), "Standard deviation of dedispersed time series")
			("nbits", value<int>()->default_value(8), "Data type of dedispersed time series")
			("savetim", "Output dedispersed data (sigproc format by default)")
			("format", value<string>()->default_value("pulsarx"), "Output format of dedispersed data [pulsarx(default),sigproc,presto]")
			("cont", "Input files are contiguous")
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

	bool contiguous = vm.count("cont");
	repeater = vm.count("repeater");

	num_threads = vm["threads"].as<unsigned int>();

	vector<double> jump = vm["jump"].as<vector<double>>();

	vector<string> fnames = vm["input"].as<vector<string>>();

	long int nfil = fnames.size();
	Filterbank *fil = new Filterbank [nfil];
	for (long int i=0; i<nfil; i++)
	{
		fil[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<nfil; i++)
	{
        fil[i].read_header();
        ntotal += fil[i].nsamples;
        MJD tstart(fil[i].tstart);
        tstarts.push_back(tstart);
        tends.push_back(tstart+fil[i].nsamples*fil[i].tsamp);
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<nfil-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*fil[idx[i]].tsamp)
		{
			if (contiguous)
			{
				cerr<<"Warning: time not contiguous"<<endl;
			}
			else
			{
				cerr<<"Error: time not contiguous"<<endl;
				exit(-1);
			}
		}
	}

	long double tstart = tstarts[idx[0]].to_day();

	string source_name = vm["source_name"].as<string>();
	string s_telescope = vm["telescope"].as<string>();
	int ibeam = vm["ibeam"].as<int>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	if (vm["source_name"].defaulted())
	{
		if (strcmp(fil[0].source_name, "") != 0)
			source_name = fil[0].source_name;
	}

	if (vm["telescope"].defaulted())
	{
		get_telescope_name(fil[0].telescope_id, s_telescope);
	}

	if (vm["ibeam"].defaulted())
	{
		if (fil[0].ibeam != 0)
			ibeam = fil[0].ibeam;
	}

	if (vm["ra"].defaulted())
	{
		if (fil[0].src_raj != 0.)
		{
			src_raj = fil[0].src_raj;
		}
	}
	if (vm["dec"].defaulted())
	{
		if (fil[0].src_dej != 0.)
		{
			src_dej = fil[0].src_dej;
		}
	}

    long int nchans = fil[0].nchans;
    double tsamp = fil[0].tsamp;
    int nifs = fil[0].nifs;

	float *buffer = new float [nchans];

	vector<SinglePulse> search1;
	parse(vm, search1);

	vector<int> tds;
	for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
	{
		tds.push_back((*sp).td*vm["td"].as<int>());
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm*td_lcm;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], fil[0].frequency_table, sizeof(double)*nchans);

	PreprocessLite prep;
	prep.td = vm["td"].as<int>();
	prep.fd = vm["fd"].as<int>();
	prep.thresig = vm["zapthre"].as<float>();
	prep.filltype = vm["fill"].as<string>();
	prep.prepare(databuf);

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	stringstream ss_ibeam;
	if (vm.count("incoherent"))
    	ss_ibeam << "ifbf" << setw(5) << setfill('0') << ibeam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << ibeam;
	string s_ibeam = ss_ibeam.str();

	long int nsearch = search1.size();
	for (long int k=0; k<nsearch; k++)
	{
		search1[k].rootname = vm["rootname"].as<string>();
		search1[k].tstart = tstart + nstart*tsamp/86400.;
		search1[k].source_name = source_name;
		search1[k].telescope = s_telescope;
		search1[k].ibeam = ibeam;
		search1[k].src_raj = src_raj;
		search1[k].src_dej = src_dej;

		search1[k].fildedisp = fil[0];
		search1[k].fildedisp.tstart = tstarts[idx[0]].to_day();
		search1[k].filltype = vm["fill"].as<string>();
		search1[k].prepare(prep);
	}

	if (verbose == 1)
	{
		cout<<"Maximum width = "<<setprecision(1)<<fixed<<search1[0].maxw*1000<<" (ms)"<<endl;
	}

	int sumif = nifs>2? 2:nifs;
	
	long int ntot = 0;
	long int ntot2 = 0;
	long int count = 0;
    long int bcnt1 = 0;
    long int bcnt2 = 0;
	for (long int idxn=0; idxn<nfil; idxn++)
	{
		long int n = idx[idxn];
        long int nseg = ceil(1.*fil[n].nsamples/NSBLK);
		long int ns_filn = 0;
		for (long int s=0; s<nseg; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/ntotal<<"%)";
			}

			fil[n].read_data(NSBLK);
#ifdef FAST
			unsigned char *pcur = (unsigned char *)(fil[n].data);
#endif
			for (long int i=0; i<NSBLK; i++)
			{
				count++;
				if (count-1<nstart or count-1>nend)
				{
					if (++ns_filn == fil[n].nsamples)
					{
						goto next;
					}
					pcur += nifs*nchans;
					continue;
				}
					
				memset(buffer, 0, sizeof(float)*nchans);
				long int m = 0;
				for (long int k=0; k<sumif; k++)
				{
					for (long int j=0; j<nchans; j++)
					{
						buffer[j] +=  pcur[m++];
					}
				}

				memcpy(&databuf.buffer[0]+bcnt1*nchans, buffer, sizeof(float)*1*nchans);
                databuf.counter ++;
				bcnt1++;
				ntot++;

				if (ntot%ndump == 0)
				{
					prep.run(databuf);
					for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
					{
						prep.isbusy = true;
						(*sp).fileid = idxn+1;
						(*sp).fname = fnames[n];
						(*sp).verbose = verbose;
						(*sp).run(prep);
					}
                    bcnt1 = 0;
				}
                
				if (++ns_filn == fil[n].nsamples)
                {
                    goto next;
                }
				pcur += nifs*nchans;
			}
		}
		next:
        fil[n].free();
	}

	std::fill(prep.buffer.begin(), prep.buffer.end(), 0.);
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

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	delete [] buffer;
	delete [] fil;

	return 0;
}
