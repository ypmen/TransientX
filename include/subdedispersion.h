/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-07-28 07:49:28
 * @modify date 2020-07-28 07:49:28
 * @desc [description]
 */

#ifndef SUBDEDISPERSION
#define SUBDEDISPERSION

#include <fstream>
#include <vector>
#include "databuffer.h"
#include "dedisperse.h"

#include "filterbank.h"

using namespace std;

#define SAVESUB 1

namespace RealTime
{
	class Subband
	{
	public:
		Subband();
		~Subband();
		void prepare();
		void run(vector<float> &data);
		void cache();
		void get_subdata(vector<float> &subdata, int idm, bool overlaped=false) const;
		void get_timdata(vector<float> &timdata, int idm, bool overlaped=false) const;
		void dumpsubdata(const string &rootname, int idm) const
		{
			ofstream outfile;
			outfile.open(rootname+".sub", ios::binary|ios::app);

			vector<float> subdata;
			get_subdata(subdata, idm);
			vector<float> subdataT(ndump*nchans, 0.);
			transpose_pad<float>(&subdataT[0], &subdata[0], nchans, ndump);

			outfile.write((char *)(&subdataT[0]), sizeof(float)*ndump*nchans);
		
			outfile.close();
		}

		void dumptimdata(const string &rootname, int idm) const
		{
			ofstream outfile;
			outfile.open(rootname+".tim", ios::binary|ios::app);

			vector<float> timdata;
			get_timdata(timdata, idm);

			outfile.write((char *)(&timdata[0]), sizeof(float)*ndump);
		
			outfile.close();
		}
	public:
		bool inplace;
		string rootname;
		int ndump;
		int noverlap;
		int nchans;
		int nsub;
		int ndm_per_sub;
		int ndm;
		long int nsamples;
		double tsamp;
		vector<double> vdm;
		vector<double> frequencies;
		vector<int> fcnt;
		vector<int> decodeidm;
		vector<int> decodeisub;
	public:
		long int counter;
		vector<int> mxdelayn;
		vector<float> buffer;
		vector<float> bufferT;
		vector<float> buffertim;
		vector<float> cachetim;
		vector<float> cachesub;
	};

	class SubbandDedispersion
	{
	public:
		SubbandDedispersion();
		SubbandDedispersion(const SubbandDedispersion &dedisp);
		SubbandDedispersion & operator=(const SubbandDedispersion &dedisp);
		~SubbandDedispersion();
		void prepare(DataBuffer<float> &databuffer);
		void run(DataBuffer<float> &databuffer, long int ns);
		void cache(){sub.cache();}
		void modifynblock();
		void makeinf(Filterbank &fil);
		void preparedump(Filterbank &fil, int nbits, const string &format);
		void rundump(float mean, float std, int nbits, const string &format);
		void get_subdata(vector<float> &subdata, int idm, bool overlaped=false) const
		{
			sub.get_subdata(subdata, idm, overlaped);
		}

		void get_timdata(vector<float> &timdata, int idm, bool overlaped=false) const
		{
			sub.get_timdata(timdata, idm, overlaped);
		}

		void dumpsubdata(const string &rootname, int idm) const
		{
			sub.dumpsubdata(rootname, idm);
		}

		void dumptimdata(const string &rootname, int idm) const
		{
			sub.dumptimdata(rootname, idm);
		}

	public:
		string rootname;
		int ndump;
		double dms;
		double ddm;
		int ndm;
		double overlap;
	public:
		float mean;
		float var;
		long int counter;
		int offset;
		int noverlap;
		int nsubband;
		int nchans;
		long int nsamples;
		double tsamp;
		vector<double> frequencies;
		vector<int> mxdelayn;
		vector<int> fmap;
		vector<int> fcnt;
		vector<double> frefsub;
		vector<float> buffer;
		vector<float> bufferT;
		vector<float> buffersub;
		vector<float> buffersubT;
		int nsub;
		long int ntot;
		Subband sub;
		std::vector<std::ofstream> outfiles;
	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
		}
	};

	struct TDMTHeader{
		uint32_t headersize; // Header size in bytes
		double tsamp; // Sampling time in seconds
		double fcentre; // Centre frequency in MHz
		double bandwidth; // Bandwidth in MHz
		double acceleration; // Acceleration in m/s/s
		uint32_t nblocks; // The number of DM-T blocks
		double dms; // The start DM in pc cm^{-3}
		double ddm; // The DM step size in pc cm^{-3}
		uint32_t ndm; // The total number of DM trials
		uint32_t blocksize; // The number of time samples per DM-T block
		uint32_t nbits; // The data encoding
		/*
		*****************************************************************
		* Dimensions of data after the header are (from outer to inner):
		* [nblocks, ndm, blocksize]
		* The value of nbits corresponds to the following encodings:
		* - 8  == int8_t
		* - 32 == float
		* Other values are currently undefined
		*****************************************************************
		*/
	   }__attribute__((packed));
}

#endif /* SUBDEDISPERSION */
