/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-03-19 19:26:11
 * @modify date 2021-03-19 19:26:11
 * @desc [description]
 */

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include "databuffer.h"

class Preprocess : public DataBuffer<float>
{
public:
	Preprocess()
	{
		td = 1;
		fd = 1;
		width = 0.;
		thresig = 3.;
	}
	~Preprocess(){}
	/* input range should be in [0, 512] */
	void prepare(DataBuffer<short> &databuffer);
	DataBuffer<float> * run(DataBuffer<short> &databuffer);
	void get_stat(DataBuffer<short> &databuffer);
	DataBuffer<float> * get(){return this;}
public:
	int td;
	int fd;
	float width;
	float thresig;
public:
	std::vector<float> chkurtosis;
	std::vector<float> chskewness;
	std::vector<float> chmean;
	std::vector<float> chstd;
};

inline void forward(const std::vector<size_t> &hist, int &median, int &res)
{
	if (res-1 < 0)
	{
		while (hist[++median] == 0);
		res = hist[median]-1;
	}
	else
	{
		res--;
	}
}

inline void backward(const std::vector<size_t> &hist, int &median, int &res)
{
	if (res+1 >= hist[median])
	{
		while (hist[--median] == 0);
		res = 0;
	}
	else
	{
		res++;
	}
}

#endif /* PREPROCESS_H */
