/*
 * downsample.h
 *
 *  Created on: May 5, 2020
 *      Author: ypmen
 */

#ifndef DOWNSAMPLE_H_
#define DOWNSAMPLE_H_

#include "databuffer.h"

class Downsample : public DataBuffer<float>
{
public:
	Downsample();
	Downsample(const Downsample &downsample);
	Downsample & operator=(const Downsample &downsample);
	Downsample(int tds, int fds);
	~Downsample();
	void prepare(DataBuffer<float> &databuffer);
	DataBuffer<float> * run(DataBuffer<float> &databuffer);
	DataBuffer<float> * get(){return this;}
public:
	int td;
	int fd;
};


#endif /* DOWNSAMPLE_H_ */