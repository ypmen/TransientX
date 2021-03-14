/*
 * cluster.cpp
 *
 *  Created on: May 1, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <vector>

#include "cluster.h"
#include "clustering.h"
#include "kdtree.h"
#include "subdedispersion.h"

using namespace std;

template <typename T>
Cluster<T>::Cluster()
{
	counter = 0;
	tsamp = 0.;
	nsamples = 0;
	dms = 0.;
	ddm = 0.;
	ndm = 0;
}

template <typename T>
Cluster<T>::Cluster(const Cluster &cluster)
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

template <typename T>
Cluster<T> & Cluster<T>::operator=(const Cluster &cluster)
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

template <typename T>
Cluster<T>::~Cluster(){}

template <typename T>
bool Cluster<T>::run(Boxcar &boxcar, float threS, double radius_smearing, int kvalue, int maxncand, int minpts, bool remove_cand_with_maxwidth)
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
				pos_smearing[0] = RealTime::SubbandDedispersion::dmdelay(dm, boxcar.fmax, boxcar.fmin)*0.5;
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

	vector<vector<long int>> state;
	//kdtree.runDBSCAN(radius_smearing*radius_smearing, kvalue);
	//kdtree.recycle(state);
	clustering<double>(state, candleft_smearing, leftS, kdtree, radius_smearing*radius_smearing, maxncand);

	vector<long int> cluster_mxS_i(kdtree.ncluster, -1);
	vector<float> cluster_maxS(kdtree.ncluster, 0.);
	candcluster.resize(kdtree.ncluster);

	for (auto i=state.begin(); i!=state.end(); ++i)
	{
		candstate.push_back(std::make_tuple(candleft_smearing[(*i)[0]][0], candleft_smearing[(*i)[0]][1], (*i)[1], (*i)[2]));
		
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

	for (long int k=0; k<kdtree.ncluster; k++)
	{
		tuple<long int, long int, int, float> cand; 
		cand = make_tuple(candleft[cluster_mxS_i[k]][0], candleft[cluster_mxS_i[k]][1], leftwn[cluster_mxS_i[k]], leftS[cluster_mxS_i[k]]);
		
		if (((!remove_cand_with_maxwidth) || (leftwn[cluster_mxS_i[k]] != boxcar.maxwn)) && candcluster[k].size()>=minpts)
		{
			candlist.push_back(cand);
		}
		else
		{
			candcluster.erase(candcluster.begin()+k);
		}
	}

	return true;
}

template <typename T>
void Cluster<T>::dumpstate2txt(const string fname)
{
	std::ofstream outfile;
	outfile.open(fname);

	for (auto i=candstate.begin(); i!=candstate.end(); ++i)
	{
		outfile<<std::get<0>(*i)<<" "<<std::get<1>(*i)<<" "<<std::get<2>(*i)<<" "<<std::get<3>(*i)<<std::endl;
	}

	outfile.close();
}

template class Cluster<double>;