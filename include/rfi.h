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
    DataBuffer<float> * zap(DataBuffer<float> &databuffer, const vector<pair<double, double>> &zaplist);
    DataBuffer<float> * zdot(DataBuffer<float> &databuffer);
    DataBuffer<float> * zero(DataBuffer<float> &databuffer);
    DataBuffer<float> * mask(DataBuffer<float> &databuffer, float threRFI2, int td, int fd);
    DataBuffer<float> * kadaneF(DataBuffer<float> &databuffer, float threRFI2, double widthlimit, int td, int fd);
    DataBuffer<float> * kadaneT(DataBuffer<float> &databuffer, float threRFI2, double bandlimit, int td, int fd);
    DataBuffer<float> * get(){return this;}
public:
    vector<int> weights;
};


#endif /* RFI_H_ */