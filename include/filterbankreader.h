/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-31 21:43:37
 * @modify date 2022-08-31 21:43:37
 * @desc [description]
 */

#ifndef FILTERBANKREADER_H
#define FILTERBANKREADER_H

#include "psrdatareader.h"
#include "filterbank.h"

class FilterbankReader : public PSRDataReader
{
public:
	FilterbankReader();
	~FilterbankReader();
	void check();
	void read_header();
	size_t read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading = false);
	MJD get_start_mjd_curfile(){return MJD(fil[idmap[ifile_cur]].tstart);}
	double get_tsamp_curfile(){return fil[idmap[ifile_cur]].tsamp;}
	size_t get_count_curfile(){return ns_filn;}
	size_t get_count(){return count;}
	size_t get_ifile(){return ifile_cur;}
	size_t get_ifile_ordered(){return idmap[ifile_cur];}
	void get_filterbank_template(Filterbank &filtem);

private:
	size_t ntot;
	size_t count;
	size_t ns_filn;
	std::vector<Filterbank> fil;
	long int ifile_cur;
	long int isubint_cur;
	long int isample_cur;
	bool update_file;
	bool update_subint;
};

#endif /* FILTERBANKREADER_H */
