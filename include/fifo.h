/*
 * fifo.h
 *
 *  Created on: Apr 25, 2020
 *      Author: ypmen
 */

#ifndef FIFO_H_
#define FIFO_H_

#include <iostream>
#include <vector>
#include <complex>

using namespace std;

/* (nsamples, nchans) */
template <typename T>
class FIFO
{
public:
	FIFO();
	FIFO(const FIFO<T> &fifo);
	FIFO<T> & operator=(const FIFO<T> &fifo);
	~FIFO();
	void resize(long int ns, long int nc);
	void dump2txt(const string &fname, long int ns);
	void dump2bin(const string &fname, long int ns);
	void dump(const string &fname, long int ns);
public:
	void write(T *data, long int ns, long int nc);
	void write_raw(T *data, long int ns, long int nc);
	void read(T *data, long int ns, long int nc);
	void readT_T(T *data, long int nc, long int ns);
	void writeT_T(T *data, long int nc, long int ns);
	void read_T(T *data, long int ns, long int nc);
	void write_T(T *data, long int ns, long int nc);
	void write_raw_T(T *data, long int ns, long int nc);
	void reset()
	{
		nw = 0;
		nr = 0;
	}
	void set_delay(long int *dy, long int nc);
	void get_mean_rms(vector<T> &mean, vector<T> &rms);
	const T * get_data_pointer()
	{
		return buffer;
	}
public:
	template <typename U>
	friend ostream& operator<< (ostream &os, const FIFO<U> &fifo);
public:
	long int ndata;
	long int *delay;
public:
	long int nsamples;
	long int nchans;
	T *buffer;
	long int nw;
	long int nr;
};


template <typename T>
ostream& operator<< (ostream &os, const FIFO<T> &fifo)
{
	os<<"buffer:"<<endl;
	for (long int j=0; j<fifo.nchans; j++)
	{
		for (long int i=0; i<fifo.nsamples; i++)
		{
			os<<fifo.buffer[i*fifo.nchans+j]<<" ";
		}
		os<<endl;
	}

	os<<"delay:"<<endl;
	for (long int j=0; j<fifo.nchans; j++)
	{
		os<<fifo.delay[j]<<" ";
	}
	os<<endl;

	return os;
}

/* (nchans, nsamples) */
template <typename T>
class FIFO_T
{
public:
	FIFO_T();
	FIFO_T(const FIFO_T<T> &fifo_t);
	FIFO_T<T> & operator=(const FIFO_T<T> &fifo_t);
	~FIFO_T();
	void resize(long int ns, long int nc);
	void dump2txt(const string &fname, long int ns);
	void dump2bin(const string &fname, long int ns);
	void dump(const string &fname, long int ns);
public:
	void write(T *data, long int nc, long int ns);
	void write_raw(T *data, long int nc, long int ns);
	void read(T *data, long int nc, long int ns);
	void reset()
	{
		nw = 0;
		nr = 0;
	}
	void set_delay(long int *dy, long int nc);
public:
	long int ndata;
	long int *delay;
public:
	long int nsamples;
	long int nchans;
	T *buffer;
	long int nw;
	long int nr;
};

template <typename T>
ostream& operator<< (ostream &os, const FIFO_T<T> &fifo)
{
	os<<"buffer:"<<endl;
	for (long int j=0; j<fifo.nchans; j++)
	{
		for (long int i=0; i<fifo.nsamples; i++)
		{
			os<<fifo.buffer[j*fifo.nsamples+i]<<" ";
		}
		os<<endl;
	}

	os<<"delay:"<<endl;
	for (long int j=0; j<fifo.nchans; j++)
	{
		os<<fifo.delay[j]<<" ";
	}
	os<<endl;

	return os;
}

#endif /* FIFO_H_ */
