/*
 * fifo.cpp
 *
 *  Created on: Apr 25, 2020
 *      Author: ypmen
 */


#include <iostream>
#include <fstream>
#include <string.h>
#include "fifo.h"
#include "dedisperse.h"

using namespace std;

/* (nsamples, nchans) */
template <typename T>
FIFO<T>::FIFO()
{
	ndata = 0;
	delay = NULL;
	nsamples = 0;
	nchans = 0;
	buffer = NULL;
	nw = 0;
	nr = 0;
}

template <typename T>
FIFO<T>::FIFO(const FIFO<T> &fifo)
{
	ndata = fifo.ndata;
	nsamples = fifo.nsamples;
	nchans = fifo.nchans;
	nw = fifo.nw;
	nr = fifo.nr;
	if (fifo.delay != NULL)
	{
		delay = new long int [nchans];
		memcpy(delay, fifo.delay, sizeof(long int)*nchans);
	}
	else
	{
		delay = NULL;
	}

	if (fifo.buffer != NULL)
	{
		buffer = new T [nsamples*nchans];
		memcpy(buffer, fifo.buffer, sizeof(T)*nsamples*nchans);
	}
	else
	{
		buffer = NULL;
	}
}

template <typename T>
FIFO<T> & FIFO<T>::operator=(const FIFO<T> &fifo)
{
	ndata = fifo.ndata;
	nsamples = fifo.nsamples;
	nchans = fifo.nchans;
	nw = fifo.nw;
	nr = fifo.nr;
	if (fifo.delay != NULL)
	{
		if (delay != NULL) delete [] delay;
		delay = new long int [nchans];
		memcpy(delay, fifo.delay, sizeof(long int)*nchans);
	}
	else
	{
		if (delay != NULL) delete [] delay;
		delay = NULL;
	}

	if (fifo.buffer != NULL)
	{
		if (buffer != NULL) delete [] buffer;
		buffer = new T [nsamples*nchans];
		memcpy(buffer, fifo.buffer, sizeof(T)*nsamples*nchans);
	}
	else
	{
		if (buffer != NULL) delete [] buffer;
		buffer = NULL;
	}

	return *this;	
}

template <typename T>
FIFO<T>::~FIFO()
{
	if (delay != NULL)
	{
		delete [] delay;
		delay = NULL;
	}

	if (buffer != NULL)
	{
		delete [] buffer;
		buffer = NULL;
	}
}

template <typename T>
void FIFO<T>::set_delay(long int *dy, long int nc)
{
	for (long int j=0; j<nc; j++)
	{
		delay[j] = dy[j];
	}
}

template <typename T>
void FIFO<T>::resize(long int ns, long int nc)
{
	if (ns != nsamples or nc != nchans)
	{
		if (buffer != NULL)
		{
			delete [] buffer;
			buffer = NULL;
		}

		if (delay != NULL)
		{
			delete [] delay;
			delay = NULL;
		}

		buffer = new T [ns*nc];
		delay = new long int [nc];
		nsamples = ns;
		nchans = nc;
	}

	memset(buffer, 0, sizeof(T)*nsamples*nchans);
	memset(delay, 0, sizeof(long int)*nchans);
	nw = 0;
	nr = 0;
	ndata = 0;
}

template <typename T>
void FIFO<T>::dump2txt(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ofstream::out);

	long int m=nr;
	for (long int i=0; i<ns; i++)
	{
		if (m >= nsamples)
		{
			m = 0;
		}
		for (long int j=0; j<nchans; j++)
		{
			outfile<<buffer[m*nchans+j]<<" ";
		}
		outfile<<endl;
		m++;
	}
	outfile.close();
}

template <typename T>
void FIFO<T>::dump2bin(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ofstream::binary);

	long int m=nr;
	for (long int i=0; i<ns; i++)
	{
		if (m >= nsamples)
		{
			m = 0;
		}
		outfile.write((char *)(buffer+m*nchans), sizeof(T)*nchans);
		m++;
	}
	outfile.close();
}

template <typename T>
void FIFO<T>::dump(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ios::binary|ios::app);

	T *data = new T [ns*nchans];
	read(data, ns, nchans);

	outfile.write((char *)data, sizeof(T)*ns*nchans);

	delete [] data;

	outfile.close();
}

template <typename T>
void FIFO<T>::get_mean_rms(vector<T> &mean, vector<T> &rms)
{
	mean.resize(nchans, 0);
	rms.resize(nchans, 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<nsamples; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			mean[j] += buffer[i*nchans+j];
			rms[j] += buffer[i*nchans+j]*buffer[i*nchans+j];
		}
	}
	
	for (long int j=0; j<nchans; j++)
	{
		mean[j] /= nsamples;
		rms[j] /= nsamples;
		rms[j] -= mean[j]*mean[j];
	}
}

/* (nsamples, nchans) */
template <typename T>
void FIFO<T>::write(T *data, long int ns, long int nc)
{
	T *pw = buffer+nw*nchans;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<ns; i++)
	{
		T *pd = data+i*nchans;
		for (long int j=0; j<nc; j++)
		{
			long int l = i-delay[j];

			if (nw+l >= 0 and nw+l < nsamples)
			{
				pw[l*nchans + j] = *pd;
			}
			else if (nw+l < 0)
			{
				pw[(l+nsamples)*nchans + j] = *pd;
			}
			else if (nw+l >= nsamples)
			{
				pw[(l-nsamples)*nchans + j] = *pd;
			}
			pd++;
		}
	}
	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template <typename T>
void FIFO<T>::write_raw(T *data, long int ns, long int nc)
{
	T *pw = buffer+nw*nchans;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<ns; i++)
	{
		T *pd = data+i*nchans;
		for (long int j=0; j<nc; j++)
		{
			if (nw+i < nsamples)
			{
				pw[i*nchans + j] = *pd;
			}
			else
			{
				pw[(i-nsamples)*nchans + j] = *pd;
			}
			pd++;
		}
	}
	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template <typename T>
void FIFO<T>::read(T *data, long int ns, long int nc)
{
	T *pr = buffer+nr*nchans;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int i=0; i<ns; i++)
	{
		T *pd = data+i*nchans;
		for (long int j=0; j<nc; j++)
		{
			if (nr+i < nsamples)
			{
				*pd = pr[i*nchans + j];
			}
			else
			{
				*pd = pr[(i-nsamples)*nchans + j];
			}
			pd++;
		}
	}
	ndata -= ns;
	nr += ns;
	nr = nr >= nsamples? nr-nsamples:nr;
}

template <typename T>
void FIFO<T>::readT_T(T *data, long int nc, long int ns)
{
	T *pr = buffer+nr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdr = pr + j*nsamples;
		for (long int i=0; i<ns; i++)
		{
			if (nr+i < nsamples)
			{
				data[j*ns+i] = pdr[i];
			}
			else
			{
				data[j*ns+i] = pdr[i-nsamples];
			}
		}
	}

	ndata -= ns;
	nr += ns;
	nr = nr >= nsamples? nr-nsamples:nr;
}

template <typename T>
void FIFO<T>::writeT_T(T *data, long int nc, long int ns)
{
	T *pw = buffer+nw;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdw = pw+j*nsamples;
		long int l = -delay[j];
		for (long int i=0; i<ns; i++)
		{
			if (nw+l >= 0 and nw+l < nsamples)
			{
				pdw[l] = data[j*ns+i];
			}
			else if (nw+l < 0)
			{
				pdw[l+nsamples] = data[j*ns+i];
			}
			else if (nw+l >= nsamples)
			{
				pdw[l-nsamples] = data[j*ns+i];
			}
			l++;
		}
	}

	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template <typename T>
void FIFO<T>::read_T(T *data, long int ns, long int nc)
{
	T *pr = buffer+nr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdr = pr + j*nsamples;
		for (long int i=0; i<ns; i++)
		{
			if (nr+i < nsamples)
			{
				data[i*nc+j] = pdr[i];
			}
			else
			{
				data[i*nc+j] = pdr[i-nsamples];
			}
		}
	}

	ndata -= ns;
	nr += ns;
	nr = nr >= nsamples? nr-nsamples:nr;
}

template <typename T>
void FIFO<T>::write_T(T *data, long int ns, long int nc)
{
	T *pw = buffer+nw;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdw = pw+j*nsamples;
		long int l = -delay[j];
		for (long int i=0; i<ns; i++)
		{
			if (nw+l >= 0 and nw+l < nsamples)
			{
				pdw[l] = data[i*nc+j];
			}
			else if (nw+l < 0)
			{
				pdw[l+nsamples] = data[i*nc+j];
			}
			else if (nw+l >= nsamples)
			{
				pdw[l-nsamples] = data[i*nc+j];
			}
			l++;
		}
	}

	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template <typename T>
void FIFO<T>::write_raw_T(T *data, long int ns, long int nc)
{
	T *pw = buffer+nw;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdw = pw+j*nsamples;
		long int l = 0;
		for (long int i=0; i<ns; i++)
		{
			if (nw+l < nsamples)
			{
				pdw[l] = data[i*nc+j];
			}
			else
			{
				pdw[l-nsamples] = data[i*nc+j];
			}
			l++;
		}
	}

	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

/* (nchans, nsamples) */

template <typename T>
FIFO_T<T>::FIFO_T()
{
	ndata = 0;
	delay = NULL;
	nsamples = 0;
	nchans = 0;
	buffer = NULL;
	nw = 0;
	nr = 0;
}

template <typename T>
FIFO_T<T>::FIFO_T(const FIFO_T<T> &fifo_t)
{
	ndata = fifo_t.ndata;
	nsamples = fifo_t.nsamples;
	nchans = fifo_t.nchans;
	nw = fifo_t.nw;
	nr = fifo_t.nr;
	if (fifo_t.delay != NULL)
	{
		delay = new long int [nchans];
		memcpy(delay, fifo_t.delay, sizeof(long int)*nchans);
	}
	else
	{
		delay = NULL;
	}

	if (fifo_t.buffer != NULL)
	{
		buffer = new T [nsamples*nchans];
		memcpy(buffer, fifo_t.buffer, sizeof(T)*nsamples*nchans);
	}
	else
	{
		buffer = NULL;
	}
}

template <typename T>
FIFO_T<T> & FIFO_T<T>::operator=(const FIFO_T<T> &fifo_t)
{
	ndata = fifo_t.ndata;
	nsamples = fifo_t.nsamples;
	nchans = fifo_t.nchans;
	nw = fifo_t.nw;
	nr = fifo_t.nr;
	if (fifo_t.delay != NULL)
	{
		if (delay != NULL) delete [] delay;
		delay = new long int [nchans];
		memcpy(delay, fifo_t.delay, sizeof(long int)*nchans);
	}
	else
	{
		if (delay != NULL) delete [] delay;
		delay = NULL;
	}

	if (fifo_t.buffer != NULL)
	{
		if (buffer != NULL) delete [] buffer;
		buffer = new T [nsamples*nchans];
		memcpy(buffer, fifo_t.buffer, sizeof(T)*nsamples*nchans);
	}
	else
	{
		if (buffer != NULL) delete [] buffer;
		buffer = NULL;
	}

	return *this;	
}

template <typename T>
FIFO_T<T>::~FIFO_T()
{
	if (delay != NULL)
	{
		delete [] delay;
		delay = NULL;
	}

	if (buffer != NULL)
	{
		delete [] buffer;
		buffer = NULL;
	}
}

template <typename T>
void FIFO_T<T>::set_delay(long int *dy, long int nc)
{
	for (long int j=0; j<nc; j++)
	{
		delay[j] = dy[j];
	}
}

template <typename T>
void FIFO_T<T>::resize(long int nc, long int ns)
{
	if (nc != nchans or ns != nsamples)
	{
		if (buffer != NULL)
		{
			delete [] buffer;
			buffer = NULL;
		}

		if (delay != NULL)
		{
			delete [] delay;
			delay = NULL;
		}

		buffer = new T [nc*ns];
		delay = new long int [nc];
		nsamples = ns;
		nchans = nc;
	}

	memset(buffer, 0, sizeof(T)*nchans*nsamples);
	memset(delay, 0, sizeof(long int)*nchans);
	nw = 0;
	nr = 0;
	ndata = 0;
}

template <typename T>
void FIFO_T<T>::dump2txt(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ofstream::out);

	for (long int j=0; j<nchans; j++)
	{
		long int m = nr;
		for (long int i=0; i<ns; i++)
		{
			if (m >= nsamples) m = 0;
			outfile<<buffer[j*nsamples+m]<<" ";
			m++;
		}
		outfile<<endl;
	}

	outfile.close();
}

template <typename T>
void FIFO_T<T>::dump2bin(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ofstream::binary);

	for (long int j=0; j<nchans; j++)
	{
		if (nr+ns > nsamples)
		{
			long int ns1 = nsamples-nr;
			long int ns2 = ns-ns1;

			outfile.write((char *)(buffer+j*nsamples+nr), sizeof(T)*ns1);
			outfile.write((char *)(buffer+j*nsamples+0), sizeof(T)*ns2);
		}
		else
		{
			outfile.write((char *)(buffer+j*nsamples+nr), sizeof(T)*ns);
		}
	}

	outfile.close();
}

template <typename T>
void FIFO_T<T>::dump(const string &fname, long int ns)
{
	ofstream outfile;
	outfile.open(fname, ios::binary|ios::app);

	T *data = new T [nchans*ns];
	read(data, nchans, ns);

	outfile.write((char *)data, sizeof(T)*nchans*ns);

	delete [] data;

	outfile.close();
}

template <typename T>
void FIFO_T<T>::read(T *data, long int nc, long int ns)
{
	T *pr = buffer+nr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdr = pr + j*nsamples;
		for (long int i=0; i<ns; i++)
		{
			if (nr+i < nsamples)
			{
				data[j*ns+i] = pdr[i];
			}
			else
			{
				data[j*ns+i] = pdr[i-nsamples];
			}
		}
	}

	ndata -= ns;
	nr += ns;
	nr = nr >= nsamples? nr-nsamples:nr;
}

template <typename T>
void FIFO_T<T>::write(T *data, long int nc, long int ns)
{
	T *pw = buffer+nw;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdw = pw+j*nsamples;
		long int l = -delay[j];
		for (long int i=0; i<ns; i++)
		{
			if (nw+l >= 0 and nw+l < nsamples)
			{
				pdw[l] = data[j*ns+i];
			}
			else if (nw+l < 0)
			{
				pdw[l+nsamples] = data[j*ns+i];
			}
			else if (nw+l >= nsamples)
			{
				pdw[l-nsamples] = data[j*ns+i];
			}
			l++;
		}
	}

	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template <typename T>
void FIFO_T<T>::write_raw(T *data, long int nc, long int ns)
{
	T *pw = buffer+nw;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nc; j++)
	{
		T *pdw = pw+j*nsamples;
		long int l = 0;
		for (long int i=0; i<ns; i++)
		{
			if (nw+l < nsamples)
			{
				pdw[l] = data[j*ns+i];
			}
			else
			{
				pdw[l-nsamples] = data[j*ns+i];
			}
			l++;
		}
	}

	ndata += ns;
	nw += ns;
	nw = nw >= nsamples? nw-nsamples:nw;
}

template class FIFO<unsigned char>;
template class FIFO<float>;
template class FIFO<double>;
template class FIFO<complex<float>>;
template class FIFO<complex<double>>;

template class FIFO_T<unsigned char>;
template class FIFO_T<float>;
template class FIFO_T<double>;
template class FIFO_T<complex<double>>;
