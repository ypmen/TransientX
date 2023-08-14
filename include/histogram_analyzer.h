#ifndef HISTOGRAM_ANALYZER_H
#define HISTOGRAM_ANALYZER_H

#include <fstream>
#include <vector>
#include <string>
#include <assert.h>

#include <databuffer.h>

// analyze output produced by histo
class HistogramAnalyzer
{
public:
	HistogramAnalyzer()
	{
		threshold = 3.;
	}
	
	HistogramAnalyzer(const std::string &fname)
	{
		threshold = 3.;

		read(fname);
	}
	
	~HistogramAnalyzer(){}

	void read(const std::string &fname)
	{
		std::ifstream infile;
		infile.open(fname, std::ios::binary);

		infile.seekg(0, std::ios::end);
		size_t fileSize = infile.tellg();
		infile.seekg(0, std::ios::beg);

		histogram.resize(fileSize / sizeof(size_t), 0);
		infile.read((char *)histogram.data(), fileSize);

		infile.close();
	}

	void get_statitstics()
	{
		size_t nchans = histogram.size() / 256;

		std::vector<size_t> counts(nchans, 0);

		mean.resize(nchans, 0.);

		for (size_t i=0; i<256; i++)
		{
			for (size_t j=0; j<nchans; j++)
			{
				mean[j] += histogram[i * nchans + j] * i;
				counts[j] += histogram[i * nchans + j];
			}
		}

		for (size_t j=0; j<nchans; j++)
		{
			mean[j] /= counts[j];
		}

		stddev.resize(nchans, 0.);
		skewness.resize(nchans, 0.);
		kurtosis.resize(nchans, 0.);
		weight.resize(nchans, 0.);

		for (size_t i=0; i<256; i++)
		{
			for (size_t j=0; j<nchans; j++)
			{
				double x = i - mean[j];
				double x2 = x * x;
				double x3 = x * x2;
				double x4 = x2 * x2;

				stddev[j] += histogram[i * nchans + j] * x2;
				skewness[j] += histogram[i * nchans + j] * x3;
				kurtosis[j] += histogram[i * nchans + j] * x4;
			}
		}

		for (size_t j=0; j<nchans; j++)
		{
			stddev[j] /= counts[j];
			stddev[j] = std::sqrt(stddev[j]);

			double stddev3 = stddev[j] * stddev[j] * stddev[j];
			double stddev4 = stddev3 * stddev[j];

			skewness[j] /= counts[j];
			skewness[j] /= stddev3;

			kurtosis[j] /= counts[j];
			kurtosis[j] /= stddev4;
			kurtosis[j] -= 3.;
		}

		/* calculate mean and std of chkurtosis and chskewness */
		std::vector<float> kurtosis_sort(kurtosis.begin(), kurtosis.end());
		std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::less<float>());
		float kurtosis_q1 = kurtosis_sort[kurtosis_sort.size()/4];
		std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::greater<float>());
		float kurtosis_q3 =kurtosis_sort[kurtosis_sort.size()/4];
		float kurtosis_R = kurtosis_q3-kurtosis_q1;

		std::vector<float> skewness_sort(skewness.begin(), skewness.end());
		std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::less<float>());
		float skewness_q1 = skewness_sort[skewness_sort.size()/4];
		std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::greater<float>());
		float skewness_q3 =skewness_sort[skewness_sort.size()/4];
		float skewness_R = skewness_q3-skewness_q1;

		size_t kill_count = 0;
		for (long int j=0; j<nchans; j++)
		{
			if (kurtosis[j]>=kurtosis_q1-threshold*kurtosis_R && \
				kurtosis[j]<=kurtosis_q3+threshold*kurtosis_R && \
				skewness[j]>=skewness_q1-threshold*skewness_R && \
				skewness[j]<=skewness_q3+threshold*skewness_R)
			{
				weight[j] = 1.;
			}
			else
			{
				kill_count++;
			}
		}

		killrate = kill_count * 1. / nchans;
	}
	
	void filter(DataBuffer<float> &databuffer)
	{
		size_t nchans = histogram.size() / 256;

		assert((nchans == databuffer.nchans) && "nchans not match");

		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans; j++)
			{
				databuffer.buffer[i*databuffer.nchans+j] = weight[j]*(databuffer.buffer[i*databuffer.nchans+j]-mean[j])/stddev[j];
			}
		}
	}

public:
	float threshold;

public:
	std::vector<double> mean;
	std::vector<double> stddev;
	std::vector<double> skewness;
	std::vector<double> kurtosis;
	std::vector<double> weight;
	double killrate;

private:
	std::vector<size_t> histogram;	
};

#endif /* HISTOGRAM_ANALYZER_H */
