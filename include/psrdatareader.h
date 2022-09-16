/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-20 22:19:22
 * @modify date 2022-08-20 22:19:22
 * @desc [description]
 */

#ifndef PSRDATAREADER_H
#define PSRDATAREADER_H

#include <vector>
#include <string>

#include "filterbank.h"
#include "databuffer.h"
#include "mjd.h"

class PSRDataReader
{
public:
	PSRDataReader()
	{
		sumif = true;
		contiguous = false;
		verbose == false;

		nsamples = 0;
		nifs = 0;
		nsblk = 0;
		nchans = 0;
		tsamp = 0.;
		skip_start = 0;
		skip_end = 0;

		apply_zero_off = true;
		apply_scloffs = false;
		apply_wts = false;

		is_end = false;
	}
	virtual ~PSRDataReader(){}
	virtual void check() = 0;
	virtual void read_header() = 0;
	virtual size_t read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading = false) = 0;
	virtual MJD get_start_mjd_curfile() = 0;
	virtual double get_tsamp_curfile() = 0;
	virtual size_t get_count_curfile() = 0;
	virtual size_t get_count() = 0;
	virtual size_t get_ifile() = 0;
	virtual size_t get_ifile_ordered() = 0;
	virtual void get_filterbank_template(Filterbank &fil) = 0;

public:
	void get_fmin_fmax(double &fmin, double &fmax)
	{
		fmin = *std::min_element(frequencies.begin(), frequencies.end());
		fmax = *std::max_element(frequencies.begin(), frequencies.end());

		double chbw = (fmax - fmin) / (nchans - 1);
		fmin -= 0.5 * chbw;
		fmax += 0.5 * chbw;
	}

	void get_fc_bw(double &fc, double &bw)
	{
		double fmin;
		double fmax;
		get_fmin_fmax(fmin, fmax);

		bw = fmax - fmin;
		fc = 0.5 * (fmin + fmax);
	}

	void get_fch1_foff(double &fch1, double &foff)
	{
		fch1 = frequencies.front();
		foff = (frequencies.back() - frequencies.front()) / (nchans - 1);
	}

	double get_chbw()
	{
		double fmin = *std::min_element(frequencies.begin(), frequencies.end());
		double fmax = *std::max_element(frequencies.begin(), frequencies.end());

		return (fmax - fmin) / (nchans - 1);
	}

public:
	std::vector<std::string> fnames;
	bool sumif;
	bool contiguous;
	bool verbose;
	size_t skip_start;
	size_t skip_end;
	bool apply_zero_off;
	bool apply_scloffs;
	bool apply_wts;

public:
	std::string telescope;
	std::string source_name;
	std::string ra;
	std::string dec;
	std::string beam;

public:
	MJD start_mjd;
	std::vector<MJD> mjd_starts;
	std::vector<MJD> mjd_ends;
	std::vector<size_t> idmap;
	size_t nsamples;
	size_t nifs;
	size_t nsblk;
	size_t nchans;
	double tsamp;
	std::vector<double> frequencies;
	bool is_end;
};

#endif /* PSRDATAREADER_H */
