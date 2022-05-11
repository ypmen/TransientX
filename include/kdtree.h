/*
 * kdtree.h
 *
 *  Created on: May 2, 2020
 *      Author: ypmen
 */

#ifndef KDTREE_H_
#define KDTREE_H_

#include <iostream>
#include <vector>
#include <vector>

#include "utils.h"

using namespace std;

template <typename T>
class KDnode
{
public:
	KDnode();
	~KDnode();
public:
	vector<T> point;
	long int index;
	long int depth;
	int clusterID;
	int flag;
	KDnode *left;
	KDnode *right;
};

template <typename T>
class KDtree
{
public:
	KDtree();
	KDtree(int k);
	~KDtree();
	KDnode<T> * newNode(const vector<T> &point, long int index, long int depth);
	void recycle(vector<vector<long int>> &state)
	{
		state.clear();
		recycleRec(root, state);
	}
private:
	KDnode<T> * buildRec(KDnode<T> *root, const vector<vector<size_t>> &idxs, const vector<vector<T>> &points, long int depth);
	KDnode<T> * insertRec(KDnode<T> *root, vector<T> &point, long int index, long int depth);
	void findNeighborsRec(KDnode<T> *root, vector<T> &point, T radius, vector<KDnode<T> *> &neighbors);
	void findDensityReachRec(vector<KDnode<T> *> &neighbors, T radius, int k, int clusterID);
	void runDBSCANRec(KDnode<T> *node, T radius, int k, int &clusterID);
	void showRec(KDnode<T> *root);
	void recycleRec(KDnode<T> *root, vector<vector<long int>> &state);

	T distence(vector<T> &point1, vector<T> &point2)
	{
		long int n = point1.size();
		T dis = 0;
		for (long int i=0; i<n; i++)
		{
			dis += (point1[i]-point2[i])*(point1[i]-point2[i]);
		}
		return dis;
	}

public:
	void build(const vector<vector<T>> & points)
	{
		vector<vector<size_t>> idxs;
		idxs.reserve(ndim);
		for (long int i=0; i<ndim; i++)
		{
			idxs.push_back(argsort(points, i));
		}
		
		root = buildRec(root, idxs, points, 0);
	}

	void insert(vector<T> &point, long int index)
	{
		root = insertRec(root, point, index, 0);
	}

	void findNeighbors(vector<T> &point, T radius, vector<KDnode<T> *> &neighbors)
	{
		neighbors.clear();
		findNeighborsRec(root, point, radius, neighbors);
	}

	void findNeighborsAdd(vector<T> &point, T radius, vector<KDnode<T> *> &neighbors)
	{
		findNeighborsRec(root, point, radius, neighbors);
	}

	void findDensityReach(KDnode<T> *node, T radius, int k, int clusterID)
	{
		vector<KDnode<T> *> neighbors;
		findNeighbors(node->point, radius, neighbors);

		//cout<<node->point[0]<<" "<<node->point[1]<<endl;

		if (node->flag == 0)
		{
			if (neighbors.size() >= k)
			{
				node->flag = 1;
				findDensityReachRec(neighbors, radius, k, clusterID);
			}
			else
			{
				node->flag = 2;
			}
		}
	}

	void runDBSCAN(T radius2, int k)
	{
		runDBSCANRec(root, radius2, k, ncluster);
	}

	void show()
	{
		showRec(root);
	}

public:
	int ndim;
	int ncluster;
	KDnode<T> *root;
};


#endif /* KDTREE_H_ */
