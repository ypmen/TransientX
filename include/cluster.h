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

template <typename T>
class Cluster
{
public:
	Cluster();
	Cluster(nlohmann::json &config);
	Cluster(const Cluster &cluster);
	Cluster & operator=(const Cluster &cluster);
	~Cluster();
	bool run(Boxcar &boxcar);
	void dumpstate2txt(const string fname);
public:
	float threS;
	double radius_smearing;
	int kvalue;
	int maxncand;
	int minpts;
	bool remove_cand_with_maxwidth;
public:
	long int counter;
	double tsamp;
	long int nsamples;
	double dms;
	double ddm;
	long int ndm;
	vector<tuple<long int, long int, int, float>> candlist;
	vector<vector<pair<long int, long int>>> candcluster;
	vector<tuple<T, T, long int, long int>> candstate;
};


#endif /* CLUSTER_H_ */
