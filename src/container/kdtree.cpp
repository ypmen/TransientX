/*
 * kdtree.cpp
 *
 *  Created on: May 2, 2020
 *      Author: ypmen
 */

#include "kdtree.h"

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

template <typename T>
KDnode<T>::KDnode()
{
	point.reserve(2);
	index = 0;
	depth = 0;
	clusterID = 0;
	flag = 0;
	left = NULL;
	right = NULL;
}

template <typename T>
KDnode<T>::~KDnode()
{
	if (left != NULL)
	{
		delete [] left;
		left = NULL;
	}

	if (right != NULL)
	{
		delete [] right;
		right = NULL;
	}
}

template <typename T>
KDtree<T>::KDtree()
{
	ndim = 0;
	ncluster = 0;
	root = NULL;
}

template <typename T>
KDtree<T>::KDtree(int k)
{
	ndim = k;
	ncluster = 0;
	root = NULL;
}

template <typename T>
KDtree<T>::~KDtree()
{
	if (root != NULL)
	{
		delete [] root;
		root = NULL;
	}
}

template <typename T>
void KDtree<T>::recycleRec(KDnode<T> *root, vector<vector<long int>> &state)
{
	if (root == NULL)
		return;

	vector<long int> stat(3);
	stat[0] = root->index;
	stat[1] = root->clusterID;
	stat[2] = root->flag;
	state.push_back(stat);

	recycleRec(root->left, state);
	recycleRec(root->right, state);
}

template <typename T>
KDnode<T> * KDtree<T>::newNode(const vector<T> &point, long int index, long int depth)
{
	KDnode<T> *node = new KDnode<T> [1];
	node->point = point;
	node->index = index;
	node->depth = depth;
	return node;
}

template <typename T>
KDnode<T> * KDtree<T>::buildRec(KDnode<T> *root, const vector<vector<size_t>> &idxs, const vector<vector<T>> &points, long int depth)
{
	if (idxs[0].empty())
		return NULL;

	int dim = depth%ndim;
	long int median = idxs[dim].size()/2;
	while (median>0 and points[idxs[dim][median-1]][dim]==points[idxs[dim][median]][dim])
	{
		median--;
	}

	root = newNode(points[idxs[dim][median]], idxs[dim][median], depth);

	vector<vector<size_t>> idxs_left(ndim);
	vector<vector<size_t>> idxs_right(ndim);

	for (long int k=0; k<ndim; k++)
	{
		for (long int i=0; i<idxs[k].size(); i++)
		{
			if (idxs[k][i] != idxs[dim][median])
			{
				if (points[idxs[k][i]][dim] < points[idxs[dim][median]][dim])
				{
					idxs_left[k].push_back(idxs[k][i]);
				}
				else
				{
					idxs_right[k].push_back(idxs[k][i]);
				}
			}
		}
	}

	root->left = buildRec(root->left, idxs_left, points, depth+1);
	root->right = buildRec(root->right, idxs_right, points, depth+1);

	return root;
}

template <typename T>
KDnode<T> * KDtree<T>::insertRec(KDnode<T> *root, vector<T> &point, long int index, long int depth)
{
	if (root == NULL)
		return newNode(point, index, depth);

	int dim = depth%ndim;

	if (point[dim] < root->point[dim])
	{
		root->left = insertRec(root->left, point, index, depth+1);
	}
	else
	{
		root->right = insertRec(root->right, point, index, depth+1);
	}

	return root;
}

template <typename T>
void KDtree<T>::findNeighborsRec(KDnode<T> *root, vector<T> &point, T radius, vector<KDnode<T> *> &neighbors)
{
	if (root == NULL)
		return;

	if (distence(root->point, point) <= radius)
	{
		neighbors.push_back(root);
	}

	bool isleft;
	int dim = root->depth%ndim;
	if (point[dim] < root->point[dim])
	{
		findNeighborsRec(root->left, point, radius, neighbors);
		isleft = true;
	}
	else
	{
		findNeighborsRec(root->right, point, radius, neighbors);
		isleft = false;
	}

	if (abs(point[dim]-root->point[dim])*abs(point[dim]-root->point[dim]) <= radius)
	{
		if (isleft)
			findNeighborsRec(root->right, point, radius, neighbors);
		else
			findNeighborsRec(root->left, point, radius, neighbors);
	}
}

// template <typename T>
// void KDtree<T>::findDensityReachRec(vector<KDnode<T> *> &neighbors, T radius, int k, int clusterID)
// {
// 	for (auto nbr=neighbors.begin(); nbr!=neighbors.end(); ++nbr)
// 	{
// 	    (*nbr)->clusterID = clusterID;
// 		if ((*nbr)->flag == 0)
// 		{
//             //cout<<(*nbr)->point[0]<<" "<<(*nbr)->point[1]<<endl;

// 			vector<KDnode<T> *> nbrs;
// 			findNeighbors((*nbr)->point, radius, nbrs);
// 			if (nbrs.size() >= k)
// 			{
// 				(*nbr)->flag = 1;
// 				findDensityReachRec(nbrs, radius, k, clusterID);
// 			}
// 			else
// 			{
// 				(*nbr)->flag = 2;
// 			}
// 		}
// 	}
// }

template <typename T>
void KDtree<T>::findDensityReachRec(vector<KDnode<T> *> &neighbors, T radius, int k, int clusterID)
{
	vector<KDnode<T> *> neighbors_temp;

	for (auto nbr=neighbors.begin(); nbr!=neighbors.end(); ++nbr)
	{
		(*nbr)->clusterID = clusterID;
		if ((*nbr)->flag == 0)
		{
			//cout<<(*nbr)->point[0]<<" "<<(*nbr)->point[1]<<endl;

			vector<KDnode<T> *> nbrs;
			findNeighbors((*nbr)->point, radius, nbrs);
			if (nbrs.size() >= k)
			{
				(*nbr)->flag = 1;
				neighbors_temp.push_back((*nbr));
			}
			else
			{
				(*nbr)->flag = 2;
			}
		}
	}

	for (auto nbr=neighbors_temp.begin(); nbr!=neighbors_temp.end(); ++nbr)
	{
		vector<KDnode<T> *> nbrs;
		findNeighbors((*nbr)->point, radius, nbrs);
		findDensityReachRec(nbrs, radius, k, clusterID);
	}
}

template <typename T>
void KDtree<T>::runDBSCANRec(KDnode<T> *node, T radius, int k, int &clusterID)
{
	if (node == NULL)
		return;

	if (node->flag == 0)
	{
		vector<KDnode<T> *> neighbors;
		findNeighbors(node->point, radius, neighbors);
		if (neighbors.size() >= k)
		{
			clusterID++;
		}
	}

	findDensityReach(node, radius, k, clusterID);
	runDBSCANRec(node->left, radius, k, clusterID);
	runDBSCANRec(node->right, radius, k, clusterID);
}

template <typename T>
void KDtree<T>::showRec(KDnode<T> *root)
{
	if (root == NULL)
		return;

	for (long int i=0; i<ndim; i++)
	{
		cout<<root->point[i]<<" ";
	}

	cout<<root->index<<" ";
	cout<<root->clusterID<< " ";
	cout<<root->flag<<endl;

	showRec(root->left);
	showRec(root->right);
}

template class KDnode<long int>;
template class KDtree<long int>;

template class KDnode<double>;
template class KDtree<double>;
