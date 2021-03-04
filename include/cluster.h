/*
 * cluster.h
 *
 *  Created on: May 1, 2020
 *      Author: ypmen
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <vector>
#include <list>
#include <tuple>
#include <utility>
#include <string>
#include "dedisperse.h"
#include "boxcar.h"

using namespace std;

class Cluster
{
public:
	Cluster();
	Cluster(const Cluster &cluster);
	Cluster & operator=(const Cluster &cluster);
	~Cluster();
	bool run(Boxcar &boxcar, float threS, double radius_smearing, int kvalue, bool remove_cand_with_maxwidth);
	void dumpstate2txt(const string fname);
public:
	long int counter;
	double tsamp;
	long int nsamples;
	double dms;
	double ddm;
	long int ndm;
	vector<tuple<long int, long int, int, float>> candlist;
	vector<vector<pair<long int, long int>>> candcluster;
	vector<vector<long int>> candstate;
};


#endif /* CLUSTER_H_ */
