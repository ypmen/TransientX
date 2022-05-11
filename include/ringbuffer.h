/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-06-16 18:29:17
 * @modify date 2021-06-16 18:29:17
 * @desc [description]
 */

#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include <vector>

class RingBuffer
{
public:
	RingBuffer()
	{
		writep = 0;
		endsample = 0;
		nsamples = 0;
		nchans = 0;
	}
	~RingBuffer(){}
	void resize(long int ns, long int nc)
	{
		data.clear();
		data.resize(ns*nc, 0.);
		nsamples = ns;
		nchans = nc;
	}
	void append(const std::vector<float> &spectra)
	{
		std::copy(spectra.begin(), spectra.end(), data.begin()+writep*nchans);
		if (++writep >= nsamples) writep = 0;
		endsample++;
	}
	void append()
	{
		if (++writep >= nsamples) writep = 0;
		endsample++;
	}
	void read(std::vector<float> &buffer, long int startsample, long int ns)
	{
		buffer.resize(ns*nchans, 0.);

		long int start = ((startsample-endsample)+writep)%nsamples;
		if (start < 0) start += nsamples;
		long int end = (start+ns)%nsamples;

		if (start < end)
		{
			std::copy(data.begin()+start*nchans, data.begin()+end*nchans, buffer.begin());
		}
		else
		{
			std::copy(data.begin()+start*nchans, data.end(), buffer.begin());
			std::copy(data.begin(), data.begin()+end*nchans, buffer.begin()+(nsamples-start)*nchans);
		}
	}
private:
	long int writep;
	long int endsample;
	long int nsamples;
	long int nchans;
	std::vector<float> data;
};

#endif /* RINGBUFFER_H */
