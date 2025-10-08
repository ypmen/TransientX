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
#include "logging.h"

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
	filltype = "mean";
	dms = 0;
	ddm = 1;
	ndm = 1000;
	overlap = 0.;
	minw = 1e-4;
	maxw = 2e-2;
	snrloss = 0.1;
	nbox = 0;
	iqr = false;
	repeater = false;
	thre = 7;
	radius_smearing = 0.003;
	kvalue = 2;
	maxncand = 100;
	minpts = 0;
	remove_cand_with_maxwidth = false;
	tstart = 0;
	ibeam = 1;
	src_raj = 0.;
	src_dej = 0.;
	
	saveimage = false;

	id = 1;
	fileid = 1;

	incoherent = false;
	verbose = false;

	outmean = 0.;
	outstd = 0.;
	outnbits = 0;
	savetim = false;
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
	zaplist_channel = sp.zaplist_channel;
	rfilist = sp.rfilist;
	filltype = sp.filltype;

	dms = sp.dms;
	ddm = sp.ddm;
	ndm = sp.ndm;
	overlap = sp.overlap;
	minw = sp.minw;
	maxw = sp.maxw;
	snrloss = sp.snrloss;
	nbox = sp.nbox;
	vwn = sp.vwn;
	iqr = sp.iqr;
	repeater = sp.repeater;

	thre = sp.thre;
	radius_smearing = sp.radius_smearing;
	kvalue = sp.kvalue;
	maxncand = sp.maxncand;
	minpts = sp.minpts;
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
	saveimage = sp.saveimage;

	incoherent = sp.incoherent;
	verbose = sp.verbose;

	outmean = sp.outmean;
	outstd = sp.outstd;
	outnbits = sp.outnbits;
	savetim = sp.savetim;
	format = sp.format;

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
	zaplist_channel = sp.zaplist_channel;
	rfilist = sp.rfilist;
	filltype = sp.filltype;

	dms = sp.dms;
	ddm = sp.ddm;
	ndm = sp.ndm;
	overlap = sp.overlap;
	minw = sp.minw;
	maxw = sp.maxw;
	snrloss = sp.snrloss;
	nbox = sp.nbox;
	vwn = sp.vwn;
	iqr = sp.iqr;
	repeater = sp.repeater;

	thre = sp.thre;
	radius_smearing = sp.radius_smearing;
	kvalue = sp.kvalue;
	maxncand = sp.maxncand;
	minpts = sp.minpts;
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
	saveimage = sp.saveimage;

	incoherent = sp.incoherent;

	verbose = sp.verbose;

	outmean = sp.outmean;
	outstd = sp.outstd;
	outnbits = sp.outnbits;
	savetim = sp.savetim;
	format = sp.format;

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
	downsample.closable = true;

	equalize.prepare(downsample);
	equalize.close();
	equalize.closable = true;

	baseline.width = bswidth;
	baseline.prepare(equalize);
	baseline.close();
	baseline.closable = true;

	rfi.filltype = filltype;
	rfi.zaplist = zaplist;
	rfi.zaplist_channel = zaplist_channel;
	rfi.rfilist = rfilist;
	rfi.thremask = threMask;
	rfi.threKadaneT = threKadaneT;
	rfi.threKadaneF = threKadaneF;
	rfi.bandlimitKT = bandlimitKT;
	rfi.widthlimit = widthlimit;
	rfi.prepare(baseline);
	rfi.close();
	rfi.closable = true;

	if (format == "presto") outnbits = 32;

	dedisp.dms = dms;
	dedisp.ddm = ddm;
	dedisp.ndm = ndm;
	dedisp.overlap = overlap;
	dedisp.ndump = rfi.nsamples;
	dedisp.rootname = rootname;
	dedisp.prepare(rfi);
	if (savetim)
		dedisp.preparedump(fildedisp, outnbits, format);
	
	boxcar.minw = minw;
	boxcar.maxw = maxw;
	boxcar.snrloss = snrloss;
	boxcar.iqr = iqr;
	boxcar.repeater = repeater;
	boxcar.prepare(dedisp);

	cluster.threS = thre;
	cluster.radius_smearing = radius_smearing;
	cluster.kvalue = kvalue;
	cluster.maxncand = maxncand;
	cluster.minpts = minpts;
	cluster.remove_cand_with_maxwidth = remove_cand_with_maxwidth;

	candplot.rootname = rootname;
	candplot.id = id;
	candplot.saveimage = saveimage;

	/** form obsinfo*/
	obsinfo["Source_name"] = source_name;

	obsinfo["Telescope"] = telescope;

	stringstream ss_tstart;
	ss_tstart << setprecision(13) << fixed << tstart;
	string s_tstart = ss_tstart.str();
	obsinfo["Tstart"] = s_tstart;

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
#if defined(HAVE_SOFA) || defined(HAVE_ERFA)
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
	DataBuffer<float> *data = downsample.run(databuffer);

	data = equalize.filter(*data);

	data = baseline.filter(*data);

	data = rfi.run(*data);
	
	if (!databuffer.isbusy) data->closable = true;
	dedisp.run(*data, data->nsamples);

	if (boxcar.run(dedisp))
	{
		if (cluster.run(boxcar))
		{
			candplot.plot(cluster, boxcar, dedisp, fileid, fname, obsinfo);
		}
	}
	
	dedisp.cache();
	if (savetim)
		dedisp.rundump(outmean, outstd, outnbits, format);
}

void parse(variables_map &vm, vector<SinglePulse> &search)
{
	SinglePulse sp;

	sp.bswidth = vm["baseline"].as<vector<float>>().back();

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
			else if (*opt == "zapchan")
			{
				std::string filename = *(opt+1);
				std::ifstream infile;
				infile.open(filename);
				if (infile.fail())
				{
					BOOST_LOG_TRIVIAL(error)<<filename<<" not exist";
					exit(-1);
				}
				int ch;
				while (infile >> ch)
				{
					sp.zaplist_channel.push_back(ch);
				}
				infile.close();

				advance(opt, 1);
			}
			else if (*opt=="zero" or *opt=="zdot")
			{
				vector<string> temp{*opt};
				sp.rfilist.push_back(temp);
			}
		}
	}

	sp.filltype = vm["fill"].as<string>();
	sp.bandlimit = vm["bandlimit"].as<double>();
	sp.bandlimitKT = vm["bandlimitKT"].as<double>();
	sp.widthlimit = vm["widthlimit"].as<double>();
	sp.threMask = vm["threMask"].as<float>();
	sp.threKadaneF = vm["threKadaneF"].as<float>();
	sp.threKadaneT = vm["threKadaneT"].as<float>();

	sp.dms = vm["dms"].as<double>();
	sp.ddm = vm["ddm"].as<double>();
	sp.ndm = vm["ndm"].as<int>();
	sp.overlap = vm["overlap"].as<double>();

	sp.minw = vm["minw"].as<float>();
	sp.snrloss = vm["snrloss"].as<float>();
	sp.maxw = vm["maxw"].as<float>();
	sp.iqr = vm.count("iqr");

	sp.thre = vm["thre"].as<float>();
	sp.radius_smearing = vm["radius"].as<double>()/1000.;
	sp.kvalue = vm["neighbors"].as<unsigned int>();
	sp.maxncand = vm["maxncand"].as<int>();
	sp.minpts = vm["minpts"].as<int>();
	sp.remove_cand_with_maxwidth = vm.count("drop");

	sp.rootname = vm["rootname"].as<string>();
	sp.saveimage = vm.count("saveimage");

	sp.incoherent = vm.count("incoherent");

	sp.outmean = vm["mean"].as<float>();
	sp.outstd = vm["std"].as<float>();
	sp.outnbits = vm["nbits"].as<int>();
	sp.savetim = vm.count("savetim");
	sp.format = vm["format"].as<string>();

	if (vm.count("ddplan"))
	{
		string filename = vm["ddplan"].as<string>();
		string line;
		ifstream ddplan(filename);
		if (ddplan.fail())
		{
			BOOST_LOG_TRIVIAL(error)<<filename<<" not exist";
			exit(-1);
		}
		int id = 0;
		while (getline(ddplan, line))
		{
			boost::trim(line);
			if(!isdigit(line[0])) continue;
			
			vector<string> parameters;
			boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);
			
			sp.td = stoi(parameters[0]);
			sp.fd = stoi(parameters[1]);
			sp.dms = stod(parameters[2]);
			sp.ddm = stod(parameters[3]);
			sp.ndm = stol(parameters[4]);
			sp.snrloss = stof(parameters[5]);
			sp.maxw = stof(parameters[6]);

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
				else if (*opt == "zapchan")
				{
					std::string filename = *(opt+1);
					std::ifstream infile;
					infile.open(filename);
					if (infile.fail())
					{
						BOOST_LOG_TRIVIAL(error)<<filename<<" not exist";
						exit(-1);
					}
					int ch;
					while (infile >> ch)
					{
						sp.zaplist_channel.push_back(ch);
					}
					infile.close();

					advance(opt, 1);
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

void parse_json(variables_map &vm, nlohmann::json &config, vector<SinglePulse> &search)
{
	SinglePulse sp;

	sp.incoherent = vm.count("incoherent");

	sp.outmean = vm["mean"].as<float>();
	sp.outstd = vm["std"].as<float>();
	sp.outnbits = vm["nbits"].as<int>();
	sp.savetim = vm.count("savetim");
	sp.format = vm["format"].as<string>();

	nlohmann::json config_ddplan = config["ddplan"];
	for (auto config=config_ddplan.begin(); config!=config_ddplan.end(); ++config)
	{
		nlohmann::json config_downsample = (*config)["downsample"];
		nlohmann::json config_baseline = (*config)["baseline"];
		nlohmann::json config_rfi = (*config)["rfi"];
		nlohmann::json config_dedisp = (*config)["subdedispersion"];
		nlohmann::json config_boxcar = (*config)["boxcar"];
		nlohmann::json config_clustering = (*config)["clustering"];
		nlohmann::json config_candplot = (*config)["candplot"];

		sp.td = config_downsample["td"];
		sp.fd = config_downsample["fd"];

		sp.bswidth = config_baseline["width"];

		sp.dms = config_dedisp["dms"];
		sp.ddm = config_dedisp["ddm"];
		sp.ndm = config_dedisp["ndm"];
		sp.overlap = config_dedisp["overlap"];

		sp.bandlimit = config_rfi["bandlimit"];
		sp.bandlimitKT = config_rfi["bandlimitKT"];
		sp.widthlimit = config_rfi["widthlimit"];
		sp.threMask = config_rfi["threMask"];
		sp.threKadaneF = config_rfi["threKadaneF"];
		sp.threKadaneT = config_rfi["threKadaneT"];

		sp.minw = config_boxcar["minw"];
		sp.snrloss = config_boxcar["snrloss"];
		sp.maxw = config_boxcar["maxw"];
		sp.iqr = config_boxcar["iqr"];

		sp.thre = config_clustering["thre"];
		sp.radius_smearing = (double)(config_clustering["radius"])/1000.;
		sp.kvalue = config_clustering["neighbors"];
		sp.maxncand = config_clustering["maxncand"];
		sp.minpts = config_clustering["minpts"];
		sp.remove_cand_with_maxwidth = config_clustering["drop"];

		sp.rootname = config_candplot["rootname"];
		sp.id = config_candplot["plan_id"];
		sp.saveimage = config_candplot["saveimage"];

		// parse zaplist
		auto config_zaplist = config_rfi["zaplist"];
		for (auto z=config_zaplist.begin(); z!=config_zaplist.end(); ++z)
		{
			sp.zaplist.push_back(std::pair<double, double>(z->front(), z->back()));
		}

		auto config_zaplist_channel = config_rfi["zaplist_channel"];
		for (auto z=config_zaplist_channel.begin(); z!=config_zaplist_channel.end(); ++z)
		{
			sp.zaplist_channel.push_back(*z);
		}

		// parse rfilist
		auto config_rfilist = config_rfi["rfilist"];
		for (auto r=config_rfilist.begin(); r!=config_rfilist.end(); ++r)
		{
			for (auto item=r->begin(); item!=r->end(); ++item)
			{
				if ((*item)=="mask" or (*item)=="kadaneF" or (*item)=="kadaneT")
				{
					sp.rfilist.push_back(std::vector<std::string>{*item, *(++item), *(++item)});
				}
				else if ((*item)=="zero" or (*item)=="zdot")
				{
					sp.rfilist.push_back(std::vector<std::string>{*item});
				}
			}
		}

		search.push_back(sp);
	}
}

