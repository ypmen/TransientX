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

class Candidate
{
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

		data.resize(npol*nchan*nbin, 0);
	}

	/* perform dedispersion */
	void dedisperse(bool coherent=false);

	/* perform dededispersion */
	void dededisperse(bool coherent=false);

	/* chop number of bins */
	void shrink_to_fit(int nwidth=20, int factor=1);

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

	/* time and channel downsample */
	void downsample(int td, int fd);

	/* normalize the data to unit variance and zero mean */
	void normalize();

	void remove_baseline(){};

	/* boxcar filter */
	void matched_filter(float snrloss=0.);

	/* save data to archive */
	void save2ar(const std::string &template_file);

	/* save data to png */
	void save2png(const std::string &rootname, float threS=7.);

	/* plot */
	void plot();

public:
	static double dmdelay(double dm, double fh, double fl)
	{
		return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
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
	float snr;
	double fl;
	double fh;
	std::string pngname;
	int planid;
	std::string rawdata;
public:
	double dm_maxsnr;
public:
	double tbin;
	int npol;
	int nchan;
	int nbin;
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
