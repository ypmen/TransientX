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
        ~SubbandDedispersion();
        void prepare(DataBuffer<float> &databuffer);
        void run(DataBuffer<float> &databuffer, long int ns);
        void cache(){sub.cache();}
        void preparedump();
        void rundump();
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
        Subband sub;
    public:
        static double dmdelay(double dm, double fh, double fl)
        {
            return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
        }
    };
}

#endif /* SUBDEDISPERSION */
