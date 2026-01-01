/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-13 21:19:42
 * @modify date 2021-04-13 21:19:42
 * @desc [description]
 */

#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <string>
#include <vector>
#include <utility>
#include <map>

#include "constants.h"

class Candidate
{
public:
	enum PolType {AABBCRCI, IQUV};
public:
	Candidate();
	~Candidate();
	void close()
	{
		data.clear();
		data.shrink_to_fit();

		dmtmap.clear();
		dmtmap.shrink_to_fit();

		dmtsnrmap.clear();
		dmtsnrmap.shrink_to_fit();

		dmtwidmap.clear();
		dmtwidmap.shrink_to_fit();

		frequencies.clear();
		frequencies.shrink_to_fit();

		profile.clear();
		profile.shrink_to_fit();
	}
	void resize(int np, int nc, int nb)
	{
		npol = np;
		nchan = nc;
		nbin = nb;

		data.resize((long int)npol*(long int)nchan*(long int)nbin, 0);
	}

	/* perform dedispersion */
	void dedisperse(bool coherent=false);

	/* perform dededispersion */
	void dededisperse(bool coherent=false);

	/* chop number of bins */
	void shrink_to_fit(int nwidth=20, bool pow2bin=false, int factor=1);

	/* convert f-t plane to dm-t plane*/
	void optimize();

	/* zap channels */
	void zap(const std::vector<std::pair<double, double>> &zaplist);

	/* zap channels by skewness and kurtosis*/
	void azap(float threshold);

	/* kadane filter */
	void kadaneF(int td, int fd, float threshold=10., int nwidth=0);

	/* clip outlier */
	void clip(int td, int fd, float threshold=10.);

	/* zero-DM matched filter */
	void zdot();

	/* calculate mean,var,skewness,kurtosis */
	void get_stats();

	/* sum polarizations */
	void sumif();

	/* tranform polarization type from AABBCRCI to IQUV */
	void AABBCRCI2IQUV();
	
	/* tranform polarization type from AABBCRCI to IQUV */
	void IQUV2AABBCRCI();

	/* time and channel downsample */
	void downsample(int td, int fd);

	/* normalize the data to unit variance and zero mean */
	void normalize();

	void remove_baseline(){};

	/* boxcar filter */
	void matched_filter(float snrloss=0.);

	/* boxcar filter in profile */
	void peak_search(float snrloss=0.);

	/* save data to archive */
	void save2ar(const std::string &template_file);

	/* save data to png */
	void save2png(const std::string &rootname, float threS=7., bool white=false);

	/* plot */
	void plot();

public:
	static double dmdelay(double dm, double fh, double fl)
	{
		return dispersion_delay(dm, fh, fl);
	}
public:
	bool captured;
	long double mjd_start;
	std::map<std::string, std::string> obsinfo;
public:
	std::string s_beam;
	int id;
	long double mjd;
	double dm;
	float width;
	float width_orig;
	float snr;
	double fl;
	double fh;
	std::string pngname;
	int planid;
	std::string rawdata;
public:
	double dm_maxsnr;
public:
	enum PolType poltype;
	double tbin;
	int npol;
	int nchan;
	long int nbin;
	std::vector<double> frequencies;
	std::vector<float> data;
public:
	std::vector<float> weights;
	std::vector<double> mean;
	std::vector<double> var;
	std::vector<float> skewness;
	std::vector<float> kurtosis;
	std::vector<float> autocorr1;
public:
	double dms;
	double ddm;
	int ndm;
	std::vector<float> dmtmap;
	std::vector<float> dmtsnrmap;
	std::vector<int> dmtwidmap;
	std::vector<float> profile;
private:
	bool isdedispersed;
	bool isnormalized;
	int maxdmid;
	int maxtid;
};

#endif /* CANDIDATE_H */
