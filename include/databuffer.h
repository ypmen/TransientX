/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-08-11 11:38:39
 * @modify date 2020-08-11 11:38:39
 * @desc [description]
 */


#ifndef DATABUFFER_H_
#define DATABUFFER_H_

#include <vector>
#include <algorithm>
#include <string.h>
#include <complex>

#ifdef __AVX2__
#include <boost/align/aligned_allocator.hpp>
#endif

using namespace std;

/* (nsamples, nchans) */
template <typename T>
class DataBuffer
{
public:
	DataBuffer();
	DataBuffer(const DataBuffer<T> &databuffer);
	DataBuffer<T> & operator=(const DataBuffer<T> &databuffer);
	DataBuffer(long int ns, int nc);
	virtual ~DataBuffer();
	virtual void prepare(DataBuffer<T> &databuffer);
	virtual DataBuffer<T> * run(DataBuffer<T> &databuffer);
	virtual DataBuffer<T> * filter(DataBuffer<T> &databuffer);
	virtual DataBuffer<T> * get(){return this;}
	void open();
	void close();
	void dump2txt(const string fname);
	void dump2bin(const string fname);
	void dump(const string fname);
	void resize(long int ns, int nc);
	void get_mean_rms(vector<T> &mean, vector<T> &var);
public:
	bool equalized;
	bool isbusy;
	bool closable;
	long int counter;
	long int nsamples;
	double tsamp;
	int nchans;
	vector<double> frequencies;
#ifdef __AVX2__
	vector<T, boost::alignment::aligned_allocator<T, 32>> buffer;
#else
	vector<T> buffer;
#endif
};

#endif /* DATABUFFER_H_ */