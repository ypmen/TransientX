/*
 * singlepulse.cpp
 *
 *  Created on: May 8, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>

#include "singlepulse.h"

using namespace std;

SinglePulse::SinglePulse()
{
    td = 1;
    fd = 1;

    bswidth = 0.;

    threMask = 7;
    bandlimit = 10;
    bandlimitKT = 10.;
    threKadaneT = 7;
    threKadaneF = 7;
    widthlimit = 10e-3;
    dms = 0;
    ddm = 1;
    ndm = 1000;
    minw = 1e-4;
    maxw = 2e-2;
    snrloss = 0.1;
    nbox = 0;
    thre = 7;
    radius_smearing = 0.003;
    kvalue = 2;
    remove_cand_with_maxwidth = false;
    tstart = 0;
    ibeam = 1;
    src_raj = 0.;
    src_dej = 0.;

    id = 1;
    fileid = 1;

    incoherent = false;
    verbose = false;
}

SinglePulse::SinglePulse(const SinglePulse &sp)
{
    downsample = sp.downsample;
    equalize = sp.equalize;
    rfi = sp.rfi;
    dedisp = sp.dedisp;
    boxcar = sp.boxcar;
    cluster = sp.cluster;
    candplot = sp.candplot;

    td = sp.td;
	fd = sp.fd;
    bswidth = sp.bswidth;

    threMask = sp.threMask;
    bandlimit = sp.bandlimit;
    bandlimitKT = sp.bandlimitKT;
    widthlimit = sp.widthlimit;
    threKadaneT = sp.threKadaneT;
    threKadaneF = sp.threKadaneF;
	zaplist = sp.zaplist;
	rfilist = sp.rfilist;

	dms = sp.dms;
	ddm = sp.ddm;
	ndm = sp.ndm;
	minw = sp.minw;
    maxw = sp.maxw;
	snrloss = sp.snrloss;
	nbox = sp.nbox;

	thre = sp.thre;
    radius_smearing = sp.radius_smearing;
    kvalue = sp.kvalue;
    remove_cand_with_maxwidth = sp.remove_cand_with_maxwidth;

	tstart = sp.tstart;
	ibeam = sp.ibeam;
    src_raj = sp.src_raj;
    src_dej = sp.src_dej;
	source_name = sp.source_name;
	rootname = sp.rootname;

    id = sp.id;

    telescope = sp.telescope;
    fileid = sp.fileid;
    fname = sp.fname;

    incoherent = sp.incoherent;
    verbose = sp.verbose;

    obsinfo = sp.obsinfo;
}

SinglePulse & SinglePulse::operator=(const SinglePulse &sp)
{
    downsample = sp.downsample;
    equalize = sp.equalize;
    rfi = sp.rfi;
    dedisp = sp.dedisp;
    boxcar = sp.boxcar;
    cluster = sp.cluster;
    candplot = sp.candplot;

    td = sp.td;
	fd = sp.fd;
    bswidth = sp.bswidth;

    threMask = sp.threMask;
    bandlimit = sp.bandlimit;
    bandlimitKT = sp.bandlimitKT;
    widthlimit = sp.widthlimit;
    threKadaneT = sp.threKadaneT;
    threKadaneF = sp.threKadaneF;
	zaplist = sp.zaplist;
	rfilist = sp.rfilist;

	dms = sp.dms;
	ddm = sp.ddm;
	ndm = sp.ndm;
	minw = sp.minw;
    maxw = sp.maxw;
	snrloss = sp.snrloss;
	nbox = sp.nbox;

	thre = sp.thre;
    radius_smearing = sp.radius_smearing;
    kvalue = sp.kvalue;
    remove_cand_with_maxwidth = sp.remove_cand_with_maxwidth;
    
	tstart = sp.tstart;
	ibeam = sp.ibeam;
    src_raj = sp.src_raj;
    src_dej = sp.src_dej;
	source_name = sp.source_name;
	rootname = sp.rootname;

    id = sp.id;

    telescope = sp.telescope;
    fileid = sp.fileid;
    fname = sp.fname;

    incoherent = sp.incoherent;

    verbose = sp.verbose;

    obsinfo = sp.obsinfo;

    return *this;
}

SinglePulse::~SinglePulse(){}

void SinglePulse::prepare(DataBuffer<float> &databuffer)
{
    downsample.td = td;
    downsample.fd = fd;
    downsample.prepare(databuffer);
    downsample.close();

    equalize.prepare(downsample);
    equalize.close();

    baseline.width = bswidth;
    baseline.prepare(equalize);
    baseline.close();

    rfi.prepare(baseline);
    rfi.close();

    dedisp.dms = dms;
    dedisp.ddm = ddm;
    dedisp.ndm = ndm;
    dedisp.ndump = rfi.nsamples;
    dedisp.rootname = rootname;
    dedisp.prepare(rfi);

    boxcar.prepare(dedisp);

    minw = minw<boxcar.tsamp ? boxcar.tsamp:minw;
	float wfactor = 1./((1.-snrloss)*(1.-snrloss));
    vwn.resize(0);
	vwn.push_back((int)round(minw/boxcar.tsamp));
	while (true)
	{
        int tmp_wn1 = vwn.back();
        int tmp_wn2 = tmp_wn1*wfactor;
		tmp_wn2 = tmp_wn2<tmp_wn1+1 ? tmp_wn1+1:tmp_wn2;
        if (tmp_wn2*boxcar.tsamp > maxw) break;
        vwn.push_back(tmp_wn2);
    }
    nbox = vwn.size();

    /** form obsinfo*/
    obsinfo["Source_name"] = source_name;

    obsinfo["Telescope"] = telescope;

    string s_ra, s_dec;
	get_s_radec(src_raj, src_dej, s_ra, s_dec);
	obsinfo["RA"] = s_ra;
	obsinfo["DEC"] = s_dec;

    stringstream ss_ibeam;
    if (incoherent)
    	ss_ibeam << "ifbf" << setw(5) << setfill('0') << ibeam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << ibeam;
	string s_ibeam = ss_ibeam.str();
	obsinfo["Beam"] = s_ibeam;

    double gl = 0., gb = 0.;
#ifdef HAVE_SOFA
	get_gl_gb(gl, gb, s_ra, s_dec);
#endif
	obsinfo["GL"] = to_string(gl);
	obsinfo["GB"] = to_string(gb);

    double ymw16_maxdm = 0.;
	ymw16_maxdm = get_maxdm_ymw16(gl, gb);

    obsinfo["MaxDM_YMW16"] = to_string(ymw16_maxdm);
}

void SinglePulse::run(DataBuffer<float> &databuffer)
{
    downsample.open();
    downsample.run(databuffer);
    databuffer.close();

    equalize.open();
    equalize.run(downsample);
    downsample.close();

    baseline.open();
    baseline.run(equalize);
    equalize.close();

    rfi.open();
    rfi.zap(baseline, zaplist);

	for (auto irfi = rfilist.begin(); irfi!=rfilist.end(); ++irfi)
	{
        if ((*irfi)[0] == "mask")
        {
            rfi.mask(rfi, threMask, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
        else if ((*irfi)[0] == "kadaneF")
        {
            rfi.kadaneF(rfi, threKadaneF*threKadaneF, widthlimit, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
        else if ((*irfi)[0] == "kadaneT")
        {
            rfi.kadaneT(rfi, threKadaneT*threKadaneT, bandlimitKT, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
		else if ((*irfi)[0] == "zdot")
        {
			rfi.zdot(rfi);
        }
		else if ((*irfi)[0] == "zero")
        {
			rfi.zero(rfi);
        }
	}
    baseline.close();

    dedisp.run(rfi, rfi.nsamples);
    rfi.close();

    if (boxcar.run(dedisp, vwn))
    {
        if (cluster.run(boxcar, thre, radius_smearing, kvalue, remove_cand_with_maxwidth))
        {
            candplot.plot(cluster, boxcar, dedisp, tstart, thre, rootname, id, fileid, fname, obsinfo);
        }
    }
    databuffer.open();
}

void parse(variables_map &vm, vector<SinglePulse> &search)
{
    SinglePulse sp;

    sp.td = vm["td"].as<int>();
	sp.fd = vm["fd"].as<int>();

    sp.bswidth = vm["baseline"].as<float>();

    vector<string> rfi_opts;
	if (vm.count("rfi"))
	{
        rfi_opts = vm["rfi"].as<vector<string>>();
        for (auto opt=rfi_opts.begin(); opt!=rfi_opts.end(); ++opt)
        {
            if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
            {
                vector<string> temp{*opt, *(opt+1), *(opt+2)};       
                sp.rfilist.push_back(temp);
                advance(opt, 2);
            }
            else if (*opt == "zap")
            {
                sp.zaplist.push_back(pair<double, double>(stod(*(opt+1)), stod(*(opt+2))));
                advance(opt, 2);
            }
            else if (*opt=="zero" or *opt=="zdot")
            {
                vector<string> temp{*opt};
                sp.rfilist.push_back(temp);
            }
        }
	}

	sp.bandlimit = vm["bandlimit"].as<double>();
    sp.bandlimitKT = vm["bandlimitKT"].as<double>();
    sp.widthlimit = vm["widthlimit"].as<double>();
    sp.threMask = vm["threMask"].as<float>();
    sp.threKadaneF = vm["threKadaneF"].as<float>();
    sp.threKadaneT = vm["threKadaneT"].as<float>();

    sp.dms = vm["dms"].as<double>();
	sp.ddm = vm["ddm"].as<double>();
	sp.ndm = vm["ndm"].as<int>();

    sp.minw = vm["minw"].as<float>();
    sp.snrloss = vm["snrloss"].as<float>();
    sp.maxw = vm["maxw"].as<float>();

    sp.thre = vm["thre"].as<float>();
    sp.radius_smearing = vm["radius"].as<double>()/1000.;
    sp.kvalue = vm["neighbors"].as<unsigned int>();
    sp.remove_cand_with_maxwidth = vm.count("drop");

	sp.rootname = vm["rootname"].as<string>();

    sp.incoherent = vm.count("incoherent");

    if (vm.count("ddplan"))
    {
        string filename = vm["ddplan"].as<string>();
        string line;
        ifstream ddplan(filename);
        int id = 0;
        while (getline(ddplan, line))
        {
            vector<string> parameters;
            boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);
            
            sp.td = stoi(parameters[0]);
            sp.fd = stoi(parameters[1]);
            sp.dms = stod(parameters[2]);
            sp.ddm = stod(parameters[3]);
            sp.ndm = stol(parameters[4]);
            sp.snrloss = stof(parameters[5]);
            sp.maxw = stof(parameters[6]);

            sp.rfilist.clear();
            for (auto opt=parameters.begin()+7; opt!=parameters.end(); ++opt)
            {
                if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
                {
                    vector<string> temp{*opt, *(opt+1), *(opt+2)};       
                    sp.rfilist.push_back(temp);
                    advance(opt, 2);
                }
                else if (*opt == "zap")
                {
                    sp.zaplist.push_back(pair<double, double>(stod(*(opt+1)), stod(*(opt+2))));
                    advance(opt, 2);
                }
                else if (*opt=="zero" or *opt=="zdot")
                {
                    vector<string> temp{*opt};
                    sp.rfilist.push_back(temp);
                }
            }

            sp.id = ++id;
            search.push_back(sp);
        }
        ddplan.close();
    }
    else
    {
        sp.id = 1;
        search.push_back(sp);
    }
}
