/*
 * candplot.h
 *
 *  Created on: May 4, 2020
 *      Author: ypmen
 */

#ifndef CANDPLOT_H_
#define CANDPLOT_H_

#include <vector>
#include <tuple>
#include <map>
#include <numeric>
#include <algorithm>

#include "cluster.h"
#include "subdedispersion.h"

class CandPlot
{
public:
	CandPlot();
	~CandPlot();
	void plot(const Cluster<double> &cluster, const Boxcar &boxcar, const RealTime::SubbandDedispersion &dedisp, double tstart, float threS, const string &rootname, int id, int fileid, std::string &fname, std::map<std::string, std::string> &obsinfo, bool saveimage=false);
private:
	vector<size_t> argsort(const vector<tuple<long int, long int, int, float>> &candlist)
	{
	vector<size_t> idx(candlist.size());
	iota(idx.begin(), idx.end(), 0);

	stable_sort(idx.begin(), idx.end(), [&candlist](size_t c1, size_t c2) {return get<3>(candlist[c1]) > get<3>(candlist[c2]);});

	return idx;
	}
public:
	static long int num_cand;
};

#endif /* CANDPLOT_H_ */