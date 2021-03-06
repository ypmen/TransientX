/*
 * dedisperse.cpp
 *
 *  Created on: Apr 24, 2020
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
#include "psrfits.h"
#include "integration.h"
#include "mjd.h"
#include "databuffer.h"
#include "singlepulse.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
unsigned int dbscan_radius;
unsigned int dbscan_k;

bool repeater;
bool dumptim=false;

int main(int argc, const char *argv[])
{
	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("duplicate", "Search on 1/2 overlap data segment")
			("drop", "Drop candidates with maximum search width")
			("iqr", "Calculate variance and mean based on IQR")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("td", value<int>()->default_value(1), "Time downsample")
			("fd", value<int>()->default_value(1), "Frequency downsample")
			("dms", value<double>()->default_value(0), "DM start")
			("ddm", value<double>()->default_value(1), "DM step")
			("ndm", value<int>()->default_value(200), "Number of DM")
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
			("baseline", value<float>()->default_value(0.1), "Remove baseline of width (s)")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
            ("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(10e-3), "Width limit of RFI kadaneF (s)")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
            ("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("source_name,s", value<string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("repeater", "Using 2D matched filter (under development)")
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
	bool duplicate = vm.count("duplicate");
	repeater = vm.count("repeater");

	num_threads = vm["threads"].as<unsigned int>();

	vector<double> jump = vm["jump"].as<vector<double>>();

	vector<string> fnames = vm["input"].as<vector<string>>();

	long int npsf = fnames.size();
	Psrfits *psf = new Psrfits [npsf];
	for (long int i=0; i<npsf; i++)
	{
		psf[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<npsf; i++)
	{
		psf[i].open();
		psf[i].primary.load(psf[i].fptr);
		psf[i].load_mode();
		psf[i].subint.load_header(psf[i].fptr);
		ntotal += psf[i].subint.nsamples;
		tstarts.push_back(psf[i].primary.start_mjd);
		tends.push_back(psf[i].primary.start_mjd+psf[i].subint.nsamples*psf[i].subint.tbin);
		psf[i].close();
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<npsf-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*psf[idx[i]].subint.tbin)
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

	psf[0].open();
	psf[0].primary.load(psf[0].fptr);
	psf[0].load_mode();
	psf[0].subint.load_header(psf[0].fptr);

	if (psf[0].mode != Integration::SEARCH)
	{
		cerr<<"Error: mode is not SEARCH"<<endl;
		exit(-1);
	}

	long double tstart = tstarts[idx[0]].to_day();

	string source_name = vm["source_name"].as<string>();
	string s_telescope = vm["telescope"].as<string>();
	int ibeam = vm["ibeam"].as<int>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	if (vm["source_name"].defaulted())
	{
		if (strcmp(psf[0].primary.src_name, "") != 0)
			source_name = psf[0].primary.src_name;
	}

	if (vm["telescope"].defaulted())
	{
		if (strcmp(psf[0].primary.telesop, "") != 0)
			s_telescope = psf[0].primary.telesop;
	}

	if (vm["ibeam"].defaulted())
	{
		if (strcmp(psf[0].primary.ibeam, "") != 0)
			ibeam = stoi(psf[0].primary.ibeam);
	}

	if (vm["ra"].defaulted())
	{
		if (strcmp(psf[0].primary.ra, "") != 0)
		{
			string ra = psf[0].primary.ra;
			ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
			src_raj = stod(ra);
		}
	}
	if (vm["dec"].defaulted())
	{
		if (strcmp(psf[0].primary.dec, "") != 0)
		{
			string dec = psf[0].primary.dec;
			dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
			src_dej = stod(dec);
		}
	}

	Integration it;
	psf[0].subint.load_integration(psf[0].fptr, 0, it);

	long int nchans = it.nchan;
    double tsamp = psf[0].subint.tbin;
    int nifs = it.npol;

	float *buffer = new float [nchans];

	vector<SinglePulse> search1;
	parse(vm, search1);
	vector<SinglePulse> search2;
	if (duplicate) parse(vm, search2);

	vector<int> tds;
	for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
	{
		tds.push_back((*sp).td);
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm*td_lcm;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], it.frequencies, sizeof(double)*nchans);

	/** duplicate */
	DataBuffer<float> databuf2;
	if (duplicate)
	{
		databuf2.resize(ndump, nchans);
		databuf2.tsamp = tsamp;
		memcpy(&databuf2.frequencies[0], it.frequencies, sizeof(double)*nchans);
	}

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	long int nsearch = search1.size();
	for (long int k=0; k<nsearch; k++)
	{
		search1[k].tstart = tstart+nstart*tsamp/86400.;
		search1[k].source_name = source_name;
		search1[k].telescope = s_telescope;
		search1[k].ibeam = ibeam;
		search1[k].src_raj = src_raj;
		search1[k].src_dej = src_dej;
		search1[k].prepare(databuf);

		if (duplicate)
		{
			search2[k].tstart = tstart+(nstart+ndump/2)*tsamp/86400.;
			search2[k].source_name = source_name;
			search2[k].telescope = s_telescope;
			search2[k].ibeam = ibeam;
			search2[k].src_raj = src_raj;
			search2[k].src_dej = src_dej;
			search2[k].prepare(databuf2);
		}
	}

	if (verbose == 1)
	{
		cout<<"Maximum width = "<<setprecision(1)<<fixed<<search1[0].maxw*1000<<" (ms)"<<endl;
	}

	psf[0].close();

	int sumif = nifs>2? 2:nifs;
	
	long int ntot = 0;
	long int ntot2 = 0;
	long int count = 0;
    long int bcnt1 = 0;
    long int bcnt2 = 0;
	for (long int idxn=0; idxn<npsf; idxn++)
	{
		long int n = idx[idxn];

		psf[n].open();
		psf[n].primary.load(psf[n].fptr);
		psf[n].load_mode();
		psf[n].subint.load_header(psf[n].fptr);

		for (long int s=0; s<psf[n].subint.nsubint; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/ntotal<<"%)";
			}

			psf[n].subint.load_integration_data(psf[n].fptr, s, it);
#ifdef FAST
			unsigned char *pcur = (unsigned char *)(it.data);
#endif
			for (long int i=0; i<it.nsblk; i++)
			{
				count++;
				if (count-1<nstart or count-1>nend)
				{
					pcur += it.npol*it.nchan;
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
                bcnt1++;
                ntot++;

				if (ntot%ndump == 0)
				{
					for (auto sp=search1.begin(); sp!=search1.end(); ++sp)
					{
						(*sp).fileid = idxn+1;
						(*sp).fname = fnames[n];
						(*sp).verbose = verbose;
						(*sp).run(databuf);
					}
                    bcnt1 = 0;
				}

				if (duplicate)
				{
					if (ntot >= ndump/2)
					{
						memcpy(&databuf2.buffer[0]+bcnt2*nchans, buffer, sizeof(float)*1*nchans);
						bcnt2++;
						ntot2++;

						if (ntot2%ndump == 0)
						{
							for (auto sp=search2.begin(); sp!=search2.end(); ++sp)
							{
								(*sp).fileid = idxn+1;
								(*sp).fname = fnames[n];
								(*sp).verbose = verbose;
								(*sp).run(databuf2);
							}
							bcnt2 = 0;
						}
					}
				}

				pcur += it.npol*it.nchan;
			}
		}
		psf[n].close();
	}

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	delete [] buffer;
	delete [] psf;

	return 0;
}
