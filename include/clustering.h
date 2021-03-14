/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-03-11 20:55:54
 * @modify date 2021-03-11 20:55:54
 * @desc [description]
 */

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <set>
#include <algorithm>
#include "kdtree.h"

/**
 * @brief 
 * 
 * @tparam T 
 * @param states id,clusterid,flag
 * @param tree 
 * @param radius2 
 * @param k 
 */
template <typename T>
void clustering(std::vector<std::vector<long int>> &states, std::vector<std::vector<T>> &points, std::vector<float> &snrs, KDtree<T> &tree, T radius2, int maxncluster=100, int k=2, int minpts=5)
{
    assert(points.size() == snrs.size());

    states.resize(points.size(), {0, 0, 0});

    std::vector<long int> idxs1;
    std::vector<long int> idxs2;
    for (long int i=0; i<points.size(); i++)
    {
        idxs1.push_back(i);
    }

    int icluster = 0;

    while ((icluster++ < maxncluster) && !idxs1.empty())
    {
        float maxsnr = -1;
        float maxsnr_idx = -1;
        for (auto idx=idxs1.begin(); idx!=idxs1.end(); ++idx)
        {
            if (snrs[*idx] > maxsnr)
            {
                maxsnr = snrs[*idx];
                maxsnr_idx = *idx;
            }
        }

        std::vector<T> point = points[maxsnr_idx];

        std::vector<KDnode<T> *> neighbors1;
        std::vector<KDnode<T> *> neighbors2;
        tree.findNeighbors(point, radius2, neighbors1);

        while (!neighbors1.empty())
        {
            for (auto nb=neighbors1.begin(); nb!=neighbors1.end(); ++nb)
            {
                long int idx = (*nb)->index;
                states[idx] = std::vector<long int>{idx, icluster, 2};
                tree.findNeighborsAdd((*nb)->point, radius2, neighbors2);
                idxs2.push_back(idx);
            }

            std::set<KDnode<T> *> neighbors_tmp(neighbors2.begin(), neighbors2.end());
            neighbors2.assign(neighbors_tmp.begin(), neighbors_tmp.end());

            neighbors1.clear();

            for (auto nb=neighbors2.begin(); nb!=neighbors2.end(); ++nb)
            {
                long int idx = (*nb)->index;
                                
                if (states[idx][2] == 0)
                {
                    neighbors1.push_back(*nb);
                }
            }

            neighbors2.clear();
        }

        std::set<long int> s(idxs2.begin(), idxs2.end());
        idxs2.assign(s.begin(), s.end());

        std::vector<long int> idxs3;
        std::set_difference(idxs1.begin(), idxs1.end(), idxs2.begin(), idxs2.end(), std::back_inserter(idxs3));

        idxs1.clear();
        idxs2.clear();
        swap(idxs1, idxs3);
    }

    tree.ncluster = icluster-1;
}

#endif /* CLUSTERING_H */
