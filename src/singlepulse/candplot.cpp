/*
 * candplot.cpp
 *
 *  Created on: May 4, 2020
 *      Author: ypmen
 */

#include <fstream>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include <string>

#include "cluster.h"
#include "boxcar.h"
#include "candplot.h"
#include "utils.h"

#include <plotx.h>

namespace plt = PlotX;

using namespace std;

long int CandPlot::num_cand = 0;

CandPlot::CandPlot(){
	rootname = "J0000-00";
	id = 0;
	saveimage = false;
}

CandPlot::CandPlot(nlohmann::json &config)
{
	rootname = config["rootname"];
	id = config["plan_id"];
	saveimage = config["saveimage"];
}

CandPlot::~CandPlot(){}

void CandPlot::plot(const Cluster<double> &cluster, const Boxcar &boxcar, const RealTime::SubbandDedispersion &dedisp, int fileid, std::string &fname, std::map<std::string, std::string> &obsinfo)
{
	vector<tuple<long int, long int, int, float>> candlist = cluster.candlist;
	vector<size_t> idx = argsort(candlist);
	long int ncand = candlist.size();

	double tstart = std::stold(obsinfo["Tstart"]);

	for (long int iidx=0; iidx<ncand; iidx++)
	{
		long int k = idx[iidx];

		long int idm = get<0>(candlist[k]);
		long int isamp = get<1>(candlist[k]);
		int wn = get<2>(candlist[k]);
		float S = get<3>(candlist[k]);
		float snr = S;
		double mjd = tstart+(-dedisp.offset+boxcar.counter-dedisp.ndump+isamp)*boxcar.tsamp/86400.;
		double bmjd = tstart+(-dedisp.offset+boxcar.counter-dedisp.ndump)*boxcar.tsamp/86400.;

		long int size = cluster.candcluster[k].size();
		long int idm_min = 100000000;
		long int idm_max = 0;
		long int isamp_min = 100000000;
		long int isamp_max = 0;
		for (long int i=0; i<size; i++)
		{
			long int idm_tmp = cluster.candcluster[k][i].first;
			idm_min = idm_tmp<idm_min ? idm_tmp:idm_min;
			idm_max = idm_tmp>idm_max ? idm_tmp:idm_max;

			long int isamp_tmp = cluster.candcluster[k][i].second;
			isamp_min = isamp_tmp<isamp_min ? isamp_tmp:isamp_min;
			isamp_max = isamp_tmp>isamp_max ? isamp_tmp:isamp_max;            
		}

		if ((idm_max-idm_min) < 20)
		{
			idm_max += 10;
			idm_min -= 10;
		}

		stringstream ss_mjd;
		ss_mjd << setprecision(10) << fixed << bmjd;
		string s_mjd = ss_mjd.str();

		stringstream ss_toa;
		ss_toa << setprecision(15) << fixed << mjd;
		string s_toa = ss_toa.str();

		stringstream ss_date;
		ss_date << setprecision(10) << fixed << tstart;
		string s_date = ss_date.str();

		string s_ibeam = obsinfo["Beam"];

		std::stringstream ss_snr;
		ss_snr<<fixed<<setprecision(1)<<snr;
		std::string s_snr = ss_snr.str();

		std::stringstream ss_width;
		ss_width<<fixed<<setprecision(1)<<wn*boxcar.tsamp*1000;
		std::string s_width = ss_width.str();

		std::stringstream ss_dm;
		ss_dm<<fixed<<setprecision(1)<<boxcar.dms+idm*boxcar.ddm;
		std::string s_dm = ss_dm.str();

		std::string s_ymw16_maxdm;
		std::stringstream ss_ymw16_maxdm;
		ss_ymw16_maxdm<<fixed<<setprecision(1)<<stod(obsinfo["MaxDM_YMW16"]);
		s_ymw16_maxdm = ss_ymw16_maxdm.str();

		double ymw16_dist = get_dist_ymw16(stod(obsinfo["GL"]), stod(obsinfo["GB"]), boxcar.dms+idm*boxcar.ddm);
		std::stringstream ss_dist;
		ss_dist<<fixed<<setprecision(1)<<ymw16_dist;
		std::string s_dist = ss_dist.str();

		stringstream ss_id;
		ss_id << setw(2) << setfill('0') << id;
		string s_id = ss_id.str();

		stringstream ss_k;
		ss_k << setw(2) << setfill('0') << iidx+1;
		string s_k = ss_k.str();

		string src_name = obsinfo["Source_name"];

		string basename = rootname + "_" + s_mjd + "_" + s_ibeam + "_" + s_k + "_" + s_id;
		string figname = basename + ".png";

		vector<float> tim, sub;
		dedisp.get_timdata(tim, idm, true);
		dedisp.get_subdata(sub, idm, true);
		long int nsamples = tim.size();
		long int nchans = sub.size()/tim.size();

		long int start_sample = isamp-15*wn < isamp_min ? isamp-15*wn:isamp_min;
		start_sample = start_sample<0 ? 0:start_sample;
		long int end_sample = isamp+15*wn > isamp_max ? isamp+15*wn:isamp_max;
		end_sample = end_sample>nsamples ? nsamples:end_sample;
		nsamples = end_sample-start_sample;

		std::vector<float> mxft;
		mxft.resize(nchans*nsamples, 0.);
		std::vector<float> vt, vf, vpf, vpfcum;
		vt.resize(nsamples, 0.);
		vf.resize(nchans, 0.);
		vpf.resize(nchans, 0.);
		vpfcum.resize(nchans, 0.);

		for (long int i=0; i<nsamples; i++)
		{
			vt[i] = (start_sample+i)*dedisp.sub.tsamp;
		}

		vector<float> vfnorm(nchans, 0.);

		float fmax = 0.;
		float fmin = 1e6;
		for (long int j=0; j<nchans; j++)
		{
			vf[j] = dedisp.sub.frequencies[j];
			vfnorm[j] = dedisp.sub.frequencies.front()+j/(nchans-1.)*(dedisp.sub.frequencies.back()-dedisp.sub.frequencies.front());
			vpf[j] = 0.;
			vpfcum[j] = 0.;
			fmax = vf[j]>fmax ? vf[j]:fmax;
			fmin = vf[j]<fmin ? vf[j]:fmin;
		}

		for (long int j=0; j<nchans; j++)
		{
			for (long int i=0; i<nsamples; i++)
			{
				mxft[j*nsamples+i] = sub[j*(dedisp.noverlap+dedisp.ndump) + start_sample+i]/dedisp.fcnt[j];
				if (i>isamp-start_sample-wn and i<isamp-start_sample+wn)
				{
					vpf[j] += mxft[j*nsamples+i];
				}
			}
		}

		float cum = 0.;
		if (vf[1] > vf[0])
		{
			for (long int j=0; j<nchans; j++)
			{
				cum += vpf[j];
				vpfcum[j] = cum;
			}
		}
		else
		{
			for (long int j=nchans-1; j>=0; j--)
			{
				cum += vpf[j];
				vpfcum[j] = cum;
			}
		}

		long int tds = wn/8;
		tds = tds<1 ? 1:tds;
		tds = tds>nsamples ? nsamples:tds;
		std::vector<float> vt_down(nsamples/tds, 0.);
		std::vector<float> mxft_down(nchans*(nsamples/tds), 0.);
		for (long int i=0; i<nsamples/tds; i++)
		{
			for (long int ii=0; ii<tds; ii++)
			{
				vt_down[i] += vt[i*tds+ii];
			}
			vt_down[i] /= tds;
		}
		for (long int j=0; j<nchans; j++)
		{
			for (long int i=0; i<nsamples/tds; i++)
			{
				for (long int ii=0; ii<tds; ii++)
				{
					mxft_down[j*(nsamples/tds)+i] += mxft[j*nsamples+i*tds+ii];
				}
				mxft_down[j*(nsamples/tds)+i] /= tds;
			}
		}

		std::vector<float> vp(nsamples, 0.);
		float vpmean = 0.;
		for (long int i=0; i<nsamples; i++)
		{
			vp[i] = tim[start_sample+i];
			vpmean += vp[i];
		}
		vpmean /= nsamples;

		float fl = 0.;
		float fh = 0.;

		long int fstart=0;
		long int fend=0;
		kadane<float>(vpf.data(), vpf.size(), &fstart, &fend);
		fl = vfnorm[fstart];
		fh = vfnorm[fend];

		float vp_down_std = 0.;
		float vp_down_mean = 0.;
		std::vector<float> vp_down(nsamples/tds, 0.);
		for (long int i=0; i<nsamples/tds; i++)
		{
			for (long int ii=0; ii<tds; ii++)
			{
				vp_down[i] += vp[i*tds+ii];
			}
			vp_down[i] /= tds;

			vp_down_mean += vp_down[i];
			vp_down_std += vp_down[i]*vp_down[i];
		}
		vp_down_mean /= nsamples/tds;
		vp_down_std /= nsamples/tds;
		vp_down_std -= vp_down_mean*vp_down_mean;
		vp_down_std = std::sqrt(vp_down_std);
		vp_down_mean = vp_down_mean == 0. ? 1. : vp_down_mean;
		vp_down_std = vp_down_std == 0. ? 1. : vp_down_std;
		for (long int i=0; i<nsamples/tds; i++)
		{
			vp_down[i] -= vp_down_mean;
			vp_down[i] /= vp_down_std;
		}

		long int ndm = boxcar.ndm;
		long int start_dm = idm-4*(idm_max-idm_min);
		long int end_dm = idm+4*(idm_max-idm_min);
		start_dm = start_dm<0 ? 0:start_dm;
		end_dm = end_dm>ndm ? ndm:end_dm;
		ndm = end_dm-start_dm;

		std::vector<float> mxdmt(ndm*nsamples, 0.);
		std::vector<float> vdm(ndm, 0.);
		std::vector<float> vSdm(ndm, 0.);
		for (long int j=0; j<ndm; j++)
		{
			vdm[j] = boxcar.dms+(start_dm+j)*boxcar.ddm;
			vSdm[j] = 0.;
		}
		for (long int j=0; j<ndm; j++)
		{
			for (long int i=0; i<nsamples; i++)
			{
				mxdmt[j*nsamples+i] = boxcar.mxS[(j+start_dm)*boxcar.nsamples+(start_sample+i)];
				vSdm[j] = mxdmt[j*nsamples+i]>vSdm[j] ? mxdmt[j*nsamples+i]:vSdm[j];
			}
		}

		int dmds = 8;
		dmds = dmds<1 ? 1:dmds;
		dmds = dmds>ndm ? ndm:dmds;

		std::vector<float> mxdmt_down((ndm/dmds)*(nsamples/tds), 0.);
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int jj=0; jj<dmds; jj++)
			{
				for (long int i=0; i<nsamples/tds; i++)
				{
					for (long int ii=0; ii<tds; ii++)
					{
						mxdmt_down[j*(nsamples/tds)+i] += mxdmt[(j*dmds+jj)*nsamples+i*tds+ii];
					}
				}
			}
		}
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int i=0; i<nsamples/tds; i++)
			{
				mxdmt_down[j*(nsamples/tds)+i] /= tds*dmds;
			}
		}

		std::vector<float> vdm_down(ndm/dmds, 0.);
		std::vector<float> vSdm_down(ndm/dmds, 0.);
		for (long int j=0; j<ndm/dmds; j++)
		{
			for (long int jj=0; jj<dmds; jj++)
			{
				vdm_down[j] += vdm[j*dmds+jj];
				vSdm_down[j] += vSdm[j*dmds+jj];
			}
			vdm_down[j] /= dmds;
			vSdm_down[j] /= dmds;
		}

		/** plot */
		plt::Figure fig(8., 1.5);

		fig.set_background_color("black");
		fig.set_default_color("white");

		float adjustx = 0., adjusty = 0.02;
		/* profile */
		plt::Axes ax_pro(0.1+adjustx, 0.75+adjustx, 0.65+adjusty, 0.80+adjusty);
		ax_pro.plot(vt_down, vp_down);
		ax_pro.annotate("S/N = "+s_snr, 0.70, 0.8);
		ax_pro.annotate("w = "+s_width+" ms", 0.70, 0.65);
		ax_pro.axvline(vt[isamp-start_sample], 0., 1., {{"color", "red"}});
		ax_pro.autoscale(true, "x", true);
		ax_pro.set_ylabel("Flux");
		ax_pro.label(true, false, false, false);
		fig.push(ax_pro);

		/* f-t */
		plt::Axes ax_ft(0.1+adjustx, 0.75+adjustx, 0.35+adjusty, 0.65+adjusty);
		ax_ft.pcolor(vt_down, vf, mxft_down, "viridis");
		ax_ft.set_ylabel("Frequency (MHz)");
		ax_ft.label(true, false, false, false);
		fig.push(ax_ft);

		/* power-f */
		plt::Axes ax_powf(0.75+adjustx, 0.95+adjustx, 0.35+adjusty, 0.65+adjusty);
		ax_powf.plot(vpfcum, vfnorm);
		ax_powf.annotate("Integral Flux", 0.1, 1.2);
		ax_powf.axhline(fl, 0., 1., {{"color", "red"}});
		ax_powf.axhline(fh, 0., 1., {{"color", "red"}});
		ax_powf.autoscale(true, "y", true);
		ax_powf.label(false, false, false, true);
		fig.push(ax_powf);

		/* dm-t */
		float xpos = vt[isamp-start_sample];
		float ypos = vdm[idm-start_dm];
		float dmpos = boxcar.dms+idm*boxcar.ddm;

		plt::Axes ax_dmt(0.1+adjustx, 0.75+adjustx, 0.05+adjusty, 0.35+adjusty);
		ax_dmt.pcolor(vt_down, vdm_down, mxdmt_down, "viridis");
		ax_dmt.circle(xpos, ypos, 5.);
		ax_dmt.set_xlabel("Time (s)");
		ax_dmt.set_ylabel("DM (cm\\u-3\\dpc)");
		fig.push(ax_dmt);

		/* dm-S */
		plt::Axes ax_dmS(0.75+adjustx, 0.95+adjustx, 0.05+adjusty, 0.35+adjusty);
		ax_dmS.plot(vSdm_down, vdm_down);
		float vSdm_down_min = *std::min_element(vSdm_down.begin(), vSdm_down.end());
		float vSdm_down_max = *std::max_element(vSdm_down.begin(), vSdm_down.end());
		float vdm_down_min = *std::min_element(vdm_down.begin(), vdm_down.end());
		ax_dmS.annotate("DM = "+s_dm+" cm\\u-3\\dpc", std::min(vSdm_down_min, cluster.threS)+0.75*(vSdm_down_max-std::min(vSdm_down_min, cluster.threS)), std::min(vdm_down_min, dmpos), {{"xycoords", "data"}, {"rotation", "270"}, {"refpos", "right"}});
		ax_dmS.axhline(dmpos, 0., 1., {{"color", "red"}});
		ax_dmS.axvline(cluster.threS, 0., 1., {{"color", "red"}});
		ax_dmS.set_xlabel("S/N");
		ax_dmS.label(false, false, true, false);
		ax_dmS.autoscale(true, "y", true);
		fig.push(ax_dmS);

		/* metadata */
		plt::Axes ax_meta(0.1+adjustx, 0.95+adjustx, 0.81+adjusty, 0.99);
		ax_meta.label(false, false, false, false);
		ax_meta.frame(false, false, false, false);
		ax_meta.minorticks_off();
		ax_meta.majorticks_off();
		std::string fontsize = "1";
		ax_meta.annotate("Telescope = " + obsinfo["Telescope"], 0.02, 0.88, {{"fontsize", fontsize}});
		ax_meta.annotate("Beam = " + obsinfo["Beam"], 0.02, 0.74, {{"fontsize", fontsize}});
		ax_meta.annotate("RA = " + obsinfo["RA"], 0.02, 0.60, {{"fontsize", fontsize}});
		ax_meta.annotate("DEC = " + obsinfo["DEC"], 0.02, 0.46, {{"fontsize", fontsize}});
		ax_meta.annotate("DM (pc/cc) = " + s_dm, 0.02, 0.32, {{"fontsize", fontsize}});
		ax_meta.annotate("Width (ms) = " + s_width, 0.02, 0.18, {{"fontsize", fontsize}});
		ax_meta.annotate("S/N = " + s_snr, 0.02, 0.04, {{"fontsize", fontsize}});
		
		ax_meta.annotate("Source name = " + obsinfo["Source_name"], 0.42, 0.88, {{"fontsize", fontsize}});
		ax_meta.annotate("Date (MJD) = " + s_toa, 0.42, 0.74, {{"fontsize", fontsize}});
		ax_meta.annotate("GL (deg) = " + obsinfo["GL"], 0.42, 0.60, {{"fontsize", fontsize}});
		ax_meta.annotate("GB (deg) = " + obsinfo["GB"], 0.42, 0.46, {{"fontsize", fontsize}});
		ax_meta.annotate("MaxDM YMW16  (pc/cc) = " + s_ymw16_maxdm, 0.42, 0.32, {{"fontsize", fontsize}});
		ax_meta.annotate("Distance YMW16 (pc) = " + s_dist, 0.42, 0.18, {{"fontsize", fontsize}});
		
		ax_meta.annotate(fname, 0.985, 0.05, {{"xycoords","figure fraction"}, {"fontsize", "0.7"}, {"rotation", "270"}});
		ax_meta.annotate("Generated by TransientX and PlotX", 0.01, 0.01, {{"xycoords","figure fraction"}, {"fontsize", "0.7"}, {"rotation", "0"}});
		fig.push(ax_meta);

		fig.save(figname+"/PNG");
		if (saveimage) fig.savepx(basename+".px");

		ofstream outfile;
		outfile.open(rootname + "_" + s_date + "_" + s_ibeam+".cands", ios_base::app); // append instead of overwrite
		
		outfile<<s_ibeam<<"\t";
		outfile<<++num_cand<<"\t";
		outfile<<setprecision(15)<<fixed<<mjd<<"\t";
		outfile<<setprecision(2)<<fixed<<boxcar.dms+idm*boxcar.ddm<<"\t";
		outfile<<setprecision(2)<<fixed<<wn*boxcar.tsamp*1000<<"\t";
		outfile<<setprecision(1)<<fixed<<snr<<"\t";
		outfile<<setprecision(0)<<fixed<<fl<<"\t";
		outfile<<setprecision(0)<<fixed<<fh<<"\t";
		outfile<<basename+".png"<<"\t";
		outfile<<s_id<<"\t";
		outfile<<fname<<endl;
	}
}
