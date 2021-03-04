/*
 * rfi.h
 *
 *  Created on: May 6, 2020
 *      Author: ypmen
 */

#ifndef RFI_H_
#define RFI_H_

#include <vector>
#include <utility>

#include "databuffer.h"

using namespace std;

class RFI : public DataBuffer<float>
{
public:
    RFI();
    RFI(const RFI &rfi);
    RFI & operator=(const RFI &rfi);
    ~RFI();
    void prepare(const DataBuffer<float> &databuffer);
    void zap(DataBuffer<float> &databuffer, const vector<pair<double, double>> &zaplist);
    void zdot(DataBuffer<float> &databuffer);
    void zero(DataBuffer<float> &databuffer);
    bool mask(DataBuffer<float> &databuffer, float threRFI2, int td, int fd);
    bool kadaneF(DataBuffer<float> &databuffer, float threRFI2, double widthlimit, int td, int fd);
    bool kadaneT(DataBuffer<float> &databuffer, float threRFI2, double bandlimit, int td, int fd);
public:
    vector<int> weights;
};


#endif /* RFI_H_ */