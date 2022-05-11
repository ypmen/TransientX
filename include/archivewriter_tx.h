/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-02 10:54:12
 * @modify date 2020-11-02 10:54:12
 * @desc [description]
 */

#ifndef ARCHIVEWRITER_TX_H
#define ARCHIVEWRITER_TX_H

#include <string>
#include <vector>

#include "psrfits.h"
#include "mjd.h"
#include "integration.h"

using namespace std;

namespace TransientX
{
	class ArchiveWriter
	{
	public:
		ArchiveWriter();
		~ArchiveWriter();
		void close();
		void prepare();
		void run(vector<float> &profile, int np, int nc, int nb, double fold_period, double tsubint, double offs_sub);
		void run(vector<float> &profiles, int ns, int np, int nc, int nb, vector<double> &fold_periods, vector<double> &tsubints, vector<double> &offs_subs);
	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
		}
	public:
		string rootname;
		string template_file;
		Integration::Mode mode;
		int ibeam;
		string src_name;
		string ra;
		string dec;
	public:
		MJD start_mjd;
		int npol;
		int nchan;
		int nbin;
		double tbin;
		double dm;
		vector<double> frequencies;
	public:
		int nsubint;
	private:
		Integration it;
		Psrfits fits;
	};
}

#endif /* ARCHIVEWRITER_H */