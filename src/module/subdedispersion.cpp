/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-07-28 09:45:43
 * @modify date 2020-07-28 09:45:43
 * @desc [description]
 */

#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h> 

#include "subdedispersion.h"
#include "logging.h"

using namespace std;
using namespace RealTime;

Subband::Subband()
{
	inplace = false;

	counter = 0;
	ndump = 0;
	noverlap = 0;
	nchans = 0;
	nsub = 0;
	ndm_per_sub = 0;
	ndm = 0;
	nsamples = 0;
	tsamp = 0.;
}

Subband::~Subband(){}

void Subband::prepare()
{
	double fmax = *max_element(frequencies.begin(), frequencies.end());

	mxdelayn.resize(nsub*nchans*ndm_per_sub, 0);
	for (long int k=0; k<nsub; k++)
	{
		for (long int j=0; j<nchans; j++)
		{
			for (long int i=0; i<ndm_per_sub; i++)
			{
				mxdelayn[k*nchans*ndm_per_sub+j*ndm_per_sub+i] = round(SubbandDedispersion::dmdelay(vdm[k*ndm_per_sub+i], fmax, frequencies[j])/tsamp);
			}
		}
	}

	buffer.resize(nsub*nsamples*nchans, 0.);
	bufferT.resize(nsub*nchans*nsamples, 0.);
	buffertim.resize(nsub*ndm_per_sub*ndump, 0.);
	cachetim.resize(nsub*ndm_per_sub*noverlap, 0.);
	cachesub.resize(nsub*nchans*(noverlap+(nsamples-ndump)), 0.);
}

void Subband::run(vector<float> &data)
{
	long int nspace = nsamples-ndump;
	for (long int k=0; k<nsub; k++)
	{
		for (long int i=0; i<ndump; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				buffer[k*nsamples*nchans+(i+nspace)*nchans+j] = data[k*ndump*nchans+i*nchans+j];
			}
		}
	}

	for (long int k=0; k<nsub; k++)
	{
		transpose_pad<float>(&bufferT[0]+k*nchans*nsamples, &buffer[0]+k*nsamples*nchans, nsamples, nchans);
	}

	fill(buffertim.begin(), buffertim.end(), 0.);
	for (long int k=0; k<nsub; k++)
	{
		for (long int j=0; j<nchans; j++)
		{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
			for (long int l=0; l<ndm_per_sub; l++)
			{
				for (long int i=0; i<ndump; i++)
				{
					buffertim[(k*ndm_per_sub+l)*ndump+i] += bufferT[(k*nchans+j)*nsamples+(i+mxdelayn[(k*nchans+j)*ndm_per_sub+l])];
				}
			}
		}
	}

	for (long int k=0; k<nsub; k++)
	{
		for (long int i=0; i<nspace; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				buffer[k*nsamples*nchans+i*nchans+j] = buffer[k*nsamples*nchans+(i+ndump)*nchans+j];
			}
		}
	}

	counter += ndump;
}

void Subband::cache()
{
	for (long int k=0; k<nsub; k++)
	{
		for (long int l=0; l<ndm_per_sub; l++)
		{
			for (long int i=0; i<noverlap; i++)
			{
				cachetim[(k*ndm_per_sub+l)*noverlap+i] = buffertim[(k*ndm_per_sub+l)*ndump+(i+ndump-noverlap)];
			}
		}
	}

	for (long int k=0; k<nsub; k++)
	{
		for (long int j=0; j<nchans; j++)
		{
			for (long int i=0; i<noverlap+nsamples-ndump; i++)
			{
				cachesub[(k*nchans+j)*(noverlap+nsamples-ndump)+i] = bufferT[(k*nchans+j)*nsamples+(i+ndump-noverlap)];
			}
		}
	}
}

void Subband::get_subdata(vector<float> &subdata, int idm, bool overlaped) const
{
	if (!overlaped)
	{
		subdata.resize(nchans*ndump, 0.);

		int isub = idm/ndm_per_sub;
		int isubdm = idm%ndm_per_sub;

		if (!inplace)
		{
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<ndump; i++)
				{
					subdata[j*ndump+i] = bufferT[(isub*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
		}
		else
		{
			int isub_real = decodeisub[isub];
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<ndump; i++)
				{
					subdata[j*ndump+i] = bufferT[(isub_real*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
		}
	}
	else
	{
		subdata.resize(nchans*(noverlap+ndump), 0.);

		int isub = idm/ndm_per_sub;
		int isubdm = idm%ndm_per_sub;

		if (!inplace)
		{
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<noverlap; i++)
				{
					subdata[j*(noverlap+ndump)+i] = cachesub[(isub*nchans+j)*(noverlap+(nsamples-ndump))+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<ndump; i++)
				{
					subdata[j*(noverlap+ndump)+(i+noverlap)] = bufferT[(isub*nchans+j)*nsamples+i+mxdelayn[(isub*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
		}
		else
		{
			int isub_real = decodeisub[isub];
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<noverlap; i++)
				{
					subdata[j*(noverlap+ndump)+i] = cachesub[(isub_real*nchans+j)*(noverlap+(nsamples-ndump))+i+mxdelayn[(isub_real*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
			for (long int j=0; j<nchans; j++)
			{
				for (long int i=0; i<ndump; i++)
				{
					subdata[j*(noverlap+ndump)+(i+noverlap)] = bufferT[(isub_real*nchans+j)*nsamples+i+mxdelayn[(isub_real*nchans+j)*ndm_per_sub+isubdm]];
				}
			}
		}
	}
}

void Subband::get_timdata(vector<float> &timdata, int idm, bool overlaped) const
{
	if (!overlaped)
	{
		timdata.resize(ndump, 0.);
		if (!inplace)
		{
			for (long int i=0; i<ndump; i++)
			{
				timdata[i] = buffertim[idm*ndump+i];
			}
		}
		else
		{
			int idm_real = decodeidm[idm];
			for (long int i=0; i<ndump; i++)
			{
				timdata[i] = buffertim[idm_real*ndump+i];
			}
		}
	}
	else
	{
		timdata.resize(noverlap+ndump, 0.);
		if (!inplace)
		{
			for (long int i=0; i<noverlap; i++)
			{
				timdata[i] = cachetim[idm*noverlap+i];
			}
			for (long int i=0; i<ndump; i++)
			{
				timdata[i+noverlap] = buffertim[idm*ndump+i];
			}
		}
		else
		{
			int idm_real = decodeidm[idm];
			for (long int i=0; i<noverlap; i++)
			{
				timdata[i] = cachetim[idm_real*noverlap+i];
			}
			for (long int i=0; i<noverlap+ndump; i++)
			{
				timdata[i+noverlap] = buffertim[idm_real*ndump+i];
			}
		}
	}
}

/** Subband dedisperison */
SubbandDedispersion::SubbandDedispersion()
{
	mean = 0.;
	var = 0.;
	counter = 0;
	offset = 0;
	noverlap = 0;
	nsubband = 0;
	ndump = 0;
	dms = 0.;
	ddm = 0.;
	ndm = 0;
	overlap = 0.;
	nchans = 0;
	nsamples = 0;
	tsamp = 0.;

	nsub = 0;
	ntot = 0;
}

SubbandDedispersion::SubbandDedispersion(const SubbandDedispersion &dedisp)
{
	rootname = dedisp.rootname;
	ndump = dedisp.ndump;
	dms = dedisp.dms;
	ddm = dedisp.ddm;
	ndm = dedisp.ndm;
	overlap = dedisp.overlap;
	mean = dedisp.mean;
	var = dedisp.var;
	counter = dedisp.counter;
	offset = dedisp.offset;
	noverlap = dedisp.noverlap;
	nsubband = dedisp.nsubband;
	nchans = dedisp.nchans;
	nsamples = dedisp.nsamples;
	tsamp = dedisp.tsamp;
	frequencies = dedisp.frequencies;
	mxdelayn = dedisp.mxdelayn;
	fmap = dedisp.fmap;
	fcnt = dedisp.fcnt;
	frefsub = dedisp.frefsub;
	buffer = dedisp.buffer;
	bufferT = dedisp.bufferT;
	buffersub = dedisp.buffersub;
	buffersubT = dedisp.buffersubT;
	nsub = dedisp.nsub;
	ntot = dedisp.ntot;
	sub = dedisp.sub;

	for (auto f=outfiles.begin(); f!=outfiles.end(); ++f)
	{
		outfiles.push_back(std::move(*f));
	}
}

SubbandDedispersion & SubbandDedispersion::operator=(const SubbandDedispersion &dedisp)
{
	rootname = dedisp.rootname;
	ndump = dedisp.ndump;
	dms = dedisp.dms;
	ddm = dedisp.ddm;
	ndm = dedisp.ndm;
	overlap = dedisp.overlap;
	mean = dedisp.mean;
	var = dedisp.var;
	counter = dedisp.counter;
	offset = dedisp.offset;
	noverlap = dedisp.noverlap;
	nsubband = dedisp.nsubband;
	nchans = dedisp.nchans;
	nsamples = dedisp.nsamples;
	tsamp = dedisp.tsamp;
	frequencies = dedisp.frequencies;
	mxdelayn = dedisp.mxdelayn;
	fmap = dedisp.fmap;
	fcnt = dedisp.fcnt;
	frefsub = dedisp.frefsub;
	buffer = dedisp.buffer;
	bufferT = dedisp.bufferT;
	buffersub = dedisp.buffersub;
	buffersubT = dedisp.buffersubT;
	nsub = dedisp.nsub;
	ntot = dedisp.ntot;
	sub = dedisp.sub;

	for (auto f=outfiles.begin(); f!=outfiles.end(); ++f)
	{
		outfiles.push_back(std::move(*f));
	}

	return *this;
}

SubbandDedispersion::~SubbandDedispersion()
{
	for (auto f=outfiles.begin(); f!=outfiles.end(); ++f)
	{
		f->close();
	}
}

void SubbandDedispersion::prepare(DataBuffer<float> &databuffer)
{
	nchans = databuffer.nchans;
	tsamp = databuffer.tsamp;
	frequencies = databuffer.frequencies;

	nsubband = round(sqrt(nchans));

	double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchans; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}

	fmap.resize(nchans, 0);
	fcnt.resize(nsubband, 0);
	frefsub.resize(nsubband, 0.);

	/** map the subband*/
	double df2 = abs(1./(fmin*fmin)-1./(fmax*fmax))/nsubband;
	double fref = frequencies[0];
	vector<double> freq(nsubband, 0.);
	long int m=0;
	for (long int i=0; i<nsubband; i++)
	{
		double fc0 = 0;
		fcnt[i] = 0;
		while(m<nchans)
		{
			double fc = frequencies[m];
			if (abs(1./(fc*fc)-1./(fref*fref))/(i+1)<=df2)
			{
				fmap[m] = i;
				fc0 += fc;
				fcnt[i]++;
				frefsub[i] = fc>frefsub[i] ? fc:frefsub[i];
				m++;
			}
			else
				break;
		}
		freq[i]=fc0/fcnt[i];
	}

	double maxsubdelayN = ceil(dmdelay(dms+ndm*ddm, fmax, fmin)/nsubband/tsamp);

	nsamples = maxsubdelayN+ndump;

	buffer.resize(nsamples*nchans, 0.);
	bufferT.resize(nchans*nsamples, 0.);
	
	/** prepare the subband */
	double ddm_sub = ddm*nsubband;
	nsub = ceil((float)ndm/nsubband);
	mxdelayn.resize(nchans*nsub, 0);
	
	noverlap = overlap*ndump;

	sub.rootname = rootname;
	sub.ndump = ndump;
	sub.noverlap = noverlap;
	sub.nchans = nsubband;
	sub.nsub = nsub;
	sub.ndm_per_sub = nsubband;
	sub.ndm = nsub*nsubband;

	sub.vdm.resize(sub.ndm, 0.);
	for (long int i=0; i<sub.ndm; i++)
	{
		sub.vdm[i] = dms + i*ddm;
	}

	double maxdelayN = ceil(dmdelay(*std::max_element(sub.vdm.begin(), sub.vdm.end()), fmax, fmin)/tsamp);
	sub.nsamples = maxdelayN+ndump;
	sub.tsamp = tsamp;

	sub.frequencies.resize(nsubband, 0.);
	sub.frequencies = frefsub;
	sub.fcnt.resize(nsubband, 0);
	sub.fcnt = fcnt;
	sub.prepare();

	/** calculate mxdelayn */
	for (long int j=0; j<nchans; j++)
	{
		for (long int k=0; k<nsub; k++)
		{
			double dm = dms + k*ddm_sub;
			mxdelayn[j*nsub+k] = round(dmdelay(dm, frefsub[fmap[j]], frequencies[j])/tsamp);
		}
	}

	buffersub.resize(nsubband*nsub*ndump, 0.);
	buffersubT.resize(nsubband*nsub*ndump, 0.);

	offset = (nsamples-ndump)+(sub.nsamples-sub.ndump)+sub.noverlap;

	std::vector<std::pair<std::string, std::string>> meta = {
		{"nsamples", std::to_string(nsamples)},
		{"ndump", std::to_string(ndump)},
		{"nchans", std::to_string(nchans)},
		{"tsamp", std::to_string(tsamp)},
		{"minimum freq", std::to_string(fmin)},
		{"maximum freq", std::to_string(fmax)},
		{"nsubband", std::to_string(nsubband)},
		{"DM start", std::to_string(dms)},
		{"DM step", std::to_string(ddm)},
		{"DM end", std::to_string(dms+ddm*ndm)},
		{"number of DM", std::to_string(ndm)},
		{"maximum delay", std::to_string(maxdelayN*tsamp)},
		{"overlap", std::to_string(overlap)}
	};
	format_logging("Subband Dedispersion Info", meta);
}

void SubbandDedispersion::run(DataBuffer<float> &databuffer, long int ns)
{
	assert(ns == ndump);

	BOOST_LOG_TRIVIAL(debug)<<"perform dedispersion on data block ("<<ndump<<", "<<nchans<<")";

	int nspace = nsamples-ns;

	for (long int i=0; i<ns; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			buffer[(i+nspace)*nchans+j] = databuffer.buffer[i*nchans+j];
		}
	}

	databuffer.isbusy = false;
	if (databuffer.closable) databuffer.close();

	transpose_pad<float>(&bufferT[0], &buffer[0], nsamples, nchans);

	fill(buffersub.begin(), buffersub.end(), 0);
	for (long int j=0; j<nchans; j++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int k=0; k<nsub; k++)
		{
			for (long int i=0; i<ndump; i++)
			{
				buffersub[fmap[j]*nsub*ndump+k*ndump+i] += bufferT[j*nsamples+i+mxdelayn[j*nsub+k]];
			}
		}
	}

	transpose_pad<float>(&buffersubT[0], &buffersub[0], nsubband, nsub*ndump);

	sub.run(buffersubT);

	for (long int i=0; i<nspace; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			buffer[i*nchans+j] = buffer[(i+ns)*nchans+j];
		}
	}

	counter += ndump;

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}

void SubbandDedispersion::preparedump(Filterbank &fil, int nbits, const string &format)
{
	double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchans; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}

	double dt = (ceil(1.*offset/ndump)*ndump-offset)*tsamp;

	outfiles.clear();
	outfiles.shrink_to_fit();

	if (format == "sigproc" or format == "presto")
	{
		struct rlimit rlim;
		int status = getrlimit(RLIMIT_NOFILE, &rlim);
		if (status)
			BOOST_LOG_TRIVIAL(warning)<<"Can't get the maximum file descriptor number";
		int rlimit_nofile = rlim.rlim_cur;
		if (ndm > rlimit_nofile)
		{
			BOOST_LOG_TRIVIAL(error)<<"Number of dedisersed files is larger than the maximum file descriptor number, please check ulimit -a";
			exit(-1);
		}

		BOOST_LOG_TRIVIAL(info)<<"create "<<ndm<<" "<<format<<" dedispersed files";
	}
	else
	{
		BOOST_LOG_TRIVIAL(info)<<"create "<<format<<" dedispersed files";
	}

	if (format == "sigproc")
	{
		for (long int k=0; k<ndm; k++)
		{
			double dm = sub.vdm[k];
			stringstream ss_dm;
			ss_dm << "DM" /*<< setw(8)*/ << setprecision(2) << fixed << setfill('0') << dm;
			string s_dm = ss_dm.str();

			std::string fname = rootname + "_" + s_dm + ".dat";
			fil.filename = fname;
			fil.tstart += dt;
			fil.refdm = dm;
			fil.tsamp = tsamp;
			fil.nchans = 1;
			fil.nbits = nbits;
			fil.nifs = 1;
			fil.data_type = 2;
			if (!fil.write_header())
				BOOST_LOG_TRIVIAL(error)<<"Error: Can not write dedisperse series header";
			fil.close();
			
			std::ofstream outfile;
			outfile.open(fname, std::ios::binary|std::ios::app);
			outfiles.push_back(std::move(outfile));
		}
	}
	else if (format == "presto")
	{
		for (long int k=0; k<ndm; k++)
		{
			double dm = sub.vdm[k];
			stringstream ss_dm;
			ss_dm << "DM" /*<< setw(8)*/ << setprecision(2) << fixed << setfill('0') << dm;
			string s_dm = ss_dm.str();
			
			std::string fname = rootname + "_" + s_dm + ".dat";
			std::ofstream outfile;
			outfile.open(fname, std::ios::binary);
			outfiles.push_back(std::move(outfile));
		}
	}
	else
	{
		TDMTHeader header;
		header.headersize = sizeof(header);
		header.tsamp = tsamp;
		header.fcentre = 0.5*(fmin+fmax);
		if (fmax != fmin)
			header.bandwidth = (fmax-fmin)/(nchans-1)*nchans;
		else
			header.bandwidth = 0.;
		header.acceleration = 0;
		header.nblocks = 0;
		header.dms = dms;
		header.ddm = ddm;
		header.ndm = ndm;
		header.blocksize = ndump;
		header.nbits = nbits;

		std::string fname = rootname+".dat";
		std::ofstream outfile;
		outfile.open(fname, std::ios::binary);
		outfiles.push_back(std::move(outfile));
		outfiles[0].write((char *)&header, sizeof(header));
		outfiles[0].close();
	}

	ntot = 0;
}

void SubbandDedispersion::modifynblock()
{
	std::string fname = rootname+".dat";
	ofstream outfile;
	outfile.open(fname, ios::in | ios::binary | ios::out);
	outfile.seekp(36);
	uint32_t nblocks = ntot/ndump;
	outfile.write((char *)&nblocks, sizeof(uint32_t));
	outfile.close();
}

void SubbandDedispersion::makeinf(Filterbank &fil)
{
	BOOST_LOG_TRIVIAL(info)<<"create presto inf files";

	double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchans; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}

	double dt = (ceil(1.*offset/ndump)*ndump-offset)*tsamp;

	for (long int k=0; k<ndm; k++)
	{
		double dm = sub.vdm[k];
		stringstream ss_dm;
		ss_dm << "DM" /*<< setw(8)*/ << setprecision(2) << fixed << setfill('0') << dm;
		string s_dm = ss_dm.str();

		std::string basename = rootname + "_" + s_dm;
		double fcentre = 0.5*(fmin+fmax);
		double bandwidth = 0.;
		if (fmax != fmin) bandwidth = (fmax-fmin)/(nchans-1)*nchans;
		double tstart = fil.tstart + dt;

		std::ofstream finf;
		finf.open(basename+".inf");
		finf<<" Data file name without suffix          =  "<<basename<<std::endl;
		finf<<" Telescope used                         =  Unknown"<<std::endl;
		finf<<" Instrument used                        =  Unknown"<<std::endl;
		finf<<" Object being observed                  =  Unknown"<<std::endl;
		finf<<" J2000 Right Ascension (hh:mm:ss.ssss)  =  00:00:00.0000"<<std::endl;
		finf<<" J2000 Declination     (dd:mm:ss.ssss)  =  00:00:00.0000"<<std::endl;
		finf<<" Data observed by                       =  Unknown"<<std::endl;
		finf<<" Epoch of observation (MJD)             =  "<<std::fixed<<std::setprecision(15)<<tstart<<std::endl;
		finf<<" Barycentered?           (1=yes, 0=no)  =  0"<<std::endl;
		finf<<" Number of bins in the time series      =  "<<ntot<<std::endl;
		finf<<" Width of each time series bin (sec)    =  "<<std::fixed<<std::setprecision(17)<<tsamp<<std::endl;
		finf<<" Any breaks in the data? (1=yes, 0=no)  =  0"<<std::endl;
		finf<<" Type of observation (EM band)          =  Radio"<<std::endl;
		finf<<" Beam diameter (arcsec)                 =  0.00"<<std::endl;
		finf<<" Dispersion measure (cm-3 pc)           =  "<<dm<<std::endl;
		finf<<" Central freq of low channel (Mhz)      =  "<<fmin<<std::endl;
		finf<<" Total bandwidth (Mhz)                  =  "<<bandwidth<<std::endl;
		finf<<" Number of channels                     =  "<<nchans<<std::endl;
		finf<<" Channel bandwidth (Mhz)                =  "<<fcentre<<std::endl;
		finf<<" Data analyzed by                       =  PulsarX"<<std::endl;
		finf.close();
	}
}

void SubbandDedispersion::rundump(float mean, float std, int nbits, const string &format)
{
	if (counter < offset+ndump) return;
	
	BOOST_LOG_TRIVIAL(debug)<<"dump dedispersed data with mean="<<mean<<" stddev="<<std;

	if (format == "sigproc" or format == "presto")
	{
		for (long int k=0; k<ndm; k++)
		{
			if (nbits == 8)
			{
				double meantim=0., stdtim=0.;
				for (long int i=0; i<sub.ndump; i++)
				{
					meantim += sub.buffertim[k*sub.ndump+i];
					stdtim += sub.buffertim[k*sub.ndump+i]*sub.buffertim[k*sub.ndump+i];
				}
				meantim /= sub.ndump;
				stdtim /= sub.ndump;
				stdtim -= meantim*meantim;
				stdtim = std::sqrt(stdtim);

				float scl = std/stdtim;
				float offs = mean-scl*meantim;

				std::vector<char> tim8bit(sub.ndump, 0);
				for (long int i=0; i<sub.ndump; i++)
				{
					float tmp = scl*sub.buffertim[k*sub.ndump+i]+offs;
					tmp = std::max(-128.f, tmp);
					tmp = std::min(127.f, tmp);
					tim8bit[i] = tmp;
				}

				outfiles[k].write((char *)(tim8bit.data()), sizeof(char)*sub.ndump);
			}
			else if (nbits == 32)
			{
				outfiles[k].write((char *)(sub.buffertim.data()+k*sub.ndump), sizeof(float)*sub.ndump);
			}
			else
			{
				outfiles[k].write((char *)(sub.buffertim.data()+k*sub.ndump), sizeof(float)*sub.ndump);
				BOOST_LOG_TRIVIAL(warning)<<"Warning: data type not supported, use float instead";
			}
		}
	}
	else
	{
		std::string fname = rootname+".dat";
		outfiles[0].open(fname, std::ios::binary|std::ios::app);

		if (nbits == 8)
		{
			std::vector<char> tim8bit(ndm*sub.ndump, 0);
			for (long int k=0; k<ndm; k++)
			{
				double meantim=0., stdtim=0.;
				for (long int i=0; i<sub.ndump; i++)
				{
					meantim += sub.buffertim[k*sub.ndump+i];
					stdtim += sub.buffertim[k*sub.ndump+i]*sub.buffertim[k*sub.ndump+i];
				}
				meantim /= sub.ndump;
				stdtim /= sub.ndump;
				stdtim -= meantim*meantim;
				stdtim = std::sqrt(stdtim);

				float scl = std/stdtim;
				float offs = mean-scl*meantim;

				for (long int i=0; i<sub.ndump; i++)
				{
					float tmp = scl*sub.buffertim[k*sub.ndump+i]+offs;
					tmp = std::max(-128.f, tmp);
					tmp = std::min(127.f, tmp);
					tim8bit[k*sub.ndump+i] = tmp;
				}
			}

			outfiles[0].write((char *)(tim8bit.data()), sizeof(unsigned char)*ndm*sub.ndump);
		}
		else if (nbits == 32)
		{
			outfiles[0].write((char *)(sub.buffertim.data()), sizeof(float)*ndm*sub.ndump);
		}
		else
		{
			outfiles[0].write((char *)(sub.buffertim.data()), sizeof(float)*ndm*sub.ndump);
			BOOST_LOG_TRIVIAL(warning)<<"Warning: data type not supported, use float instead";
		}

		outfiles[0].close();
	}

	ntot += ndump;
}
