/*
 * cluster.cpp
 *
 *  Created on: May 1, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <vector>

#include "cluster.h"
#include "kdtree.h"
#include "subdedispersion.h"

using namespace std;

Cluster::Cluster()
{
	counter = 0;
	tsamp = 0.;
	nsamples = 0;
	dms = 0.;
	ddm = 0.;
	ndm = 0;
}

Cluster::Cluster(const Cluster &cluster)
{
	counter = cluster.counter;
	tsamp = cluster.tsamp;
	nsamples = cluster.nsamples;
	dms = cluster.dms;
	ddm = cluster.ddm;
	ndm = cluster.ndm;
	candlist = cluster.candlist;
	candcluster = cluster.candcluster;
	candstate = cluster.candstate;
}

Cluster & Cluster::operator=(const Cluster &cluster)
{
	counter = cluster.counter;
	tsamp = cluster.tsamp;
	nsamples = cluster.nsamples;
	dms = cluster.dms;
	ddm = cluster.ddm;
	ndm = cluster.ndm;
	candlist = cluster.candlist;
	candcluster = cluster.candcluster;
	candstate = cluster.candstate;

	return *this;
}

Cluster::~Cluster(){}

bool Cluster::run(Boxcar &boxcar, float threS, double radius_smearing, int kvalue, bool remove_cand_with_maxwidth)
{
	counter = boxcar.counter;
	if (counter <= 0) return false;

	tsamp = boxcar.tsamp;

	candlist.clear();
	candcluster.clear();
	candstate.clear();

	vector<vector<long int>> candleft;
	vector<vector<double>> candleft_smearing;
	vector<int> leftwn;
	vector<float> leftS;

    candleft.reserve(8192);
	candleft_smearing.reserve(8192);
    leftwn.reserve(8192);
    leftS.reserve(8192);

	//remove all points smaller than the threshold
	vector<long int> pos;
	pos.resize(2);
	vector<double> pos_smearing;
	pos_smearing.resize(2);
	float *pS = boxcar.mxS;
	int *pwn = boxcar.mxwn;
	for (long int j=0; j<boxcar.ndm; j++)
	{
		for (long int i=0; i<boxcar.nsamples; i++)
		{
			if (*pS > threS)
			{
				pos[0] = j;
				pos[1] = i;
				candleft.push_back(pos);

				double dm = boxcar.dms+j*boxcar.ddm;
				pos_smearing[0] = RealTime::SubbandDedispersion::dmdelay(dm, boxcar.fmax, boxcar.fmin);
				pos_smearing[1] = i*tsamp;
				candleft_smearing.push_back(pos_smearing);

				leftwn.push_back(*pwn);
				leftS.push_back(*pS);
			}
			pS++;
			pwn++;
		}
	}

	KDtree<double> kdtree(2);
	kdtree.build(candleft_smearing);

	kdtree.runDBSCAN(radius_smearing*radius_smearing, kvalue);

	vector<vector<long int>> state;
	kdtree.recycle(state);

	vector<long int> cluster_mxS_i(kdtree.ncluster, -1);
	vector<float> cluster_maxS(kdtree.ncluster, 0.);
	candcluster.resize(kdtree.ncluster);
	vector<long int> cstate(4);

	for (auto i=state.begin(); i!=state.end(); ++i)
	{
		cstate[0] = candleft[(*i)[0]][0];
		cstate[1] = candleft[(*i)[0]][1];
		cstate[2] = (*i)[1];
		cstate[3] = (*i)[2];
		candstate.push_back(cstate);
		if ((*i)[1] !=0 )
		{
			pair<long int, long int> cand(candleft[(*i)[0]][0], candleft[(*i)[0]][1]);
			candcluster[(*i)[1]-1].push_back(cand);

			if (leftS[(*i)[0]] > cluster_maxS[(*i)[1]-1])
			{
				cluster_maxS[(*i)[1]-1] = leftS[(*i)[0]];
				cluster_mxS_i[(*i)[1]-1] = (*i)[0];
			}
		}
	}

	for (auto i=cluster_mxS_i.begin(); i!=cluster_mxS_i.end(); ++i)
	{
		tuple<long int, long int, int, float> cand; 
		cand = make_tuple(candleft[*i][0], candleft[*i][1], leftwn[*i], leftS[*i]);
		
		if ((!remove_cand_with_maxwidth) || (leftwn[*i] != boxcar.maxwn))
		{
			candlist.push_back(cand);
		}
	}

	return true;
}

void Cluster::dumpstate2txt(const string fname)
{
	FILE *fptr = fopen(fname.c_str(), "w");

	for (auto i=candstate.begin(); i!=candstate.end(); ++i)
	{
		for (auto j=(*i).begin(); j!=(*i).end(); j++)
		{
			fprintf(fptr,"%ld ", *j);
		}
		fprintf(fptr,"\n");
	}

	fclose(fptr);
}
