/*
 * singlepulse.h
 *
 *  Created on: May 8, 2020
 *      Author: ypmen
 */

#ifndef SINGLEPULSE_H_
#define SINGLEPULSEE_H_

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 

#include <string.h>

#include "databuffer.h"
#include "downsample.h"
#include "equalize.h"
#include "rfi.h"
#include "subdedispersion.h"
#include "boxcar.h"
#include "cluster.h"
#include "candplot.h"
#include "baseline.h"

using namespace std;
using namespace boost::program_options;

class SinglePulse
{
public:
    SinglePulse();
    SinglePulse(const SinglePulse &sp);
    SinglePulse & operator=(const SinglePulse &sp);
    ~SinglePulse();
    void prepare(DataBuffer<float> &databuffer);
    void run(DataBuffer<float> &databuffer);
public:
    //components
    Downsample downsample;
    Equalize equalize;
    BaseLine baseline;
    RFI rfi;
    RealTime::SubbandDedispersion dedisp;
    
    Boxcar boxcar;
    Cluster cluster;
    CandPlot candplot;

    //downsample
    int td;
    int fd;

    //baseline
    float bswidth;

    //rfi
    vector<pair<double, double>> zaplist;
    vector<vector<string>> rfilist;
    double bandlimit;
    double bandlimitKT;
    double widthlimit;
    float threKadaneT;
    float threKadaneF;
    float threMask;

    //dedispere
    double dms;
    double ddm;
    long int ndm;
    double overlap;

    //boxcar
    float minw;
    float maxw;
    float snrloss;
    int nbox;
    bool iqr;
    vector<int> vwn;

    //cluster
    float thre;
    double radius_smearing;
	int kvalue;
    bool remove_cand_with_maxwidth;

    //candplot
    long double tstart;
    int ibeam;
    double src_raj;
    double src_dej;
    string source_name;
    string rootname;
    int id;
    string telescope;

    int fileid;
    string fname;
    bool incoherent;

    bool verbose;

    std::map<std::string, std::string> obsinfo;
};

void parse(variables_map &vm, vector<SinglePulse> &search);

#endif /* SINGLEPULSE_H_ */