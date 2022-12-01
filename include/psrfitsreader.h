/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-20 22:11:40
 * @modify date 2022-08-20 22:11:40
 * @desc [description]
 */

#ifndef PSRFITSREADER_H
#define PSRFITSREADER_H

#include "psrdatareader.h"
#include "psrfits.h"

class PsrfitsReader : public PSRDataReader
{
public:
	PsrfitsReader();
	~PsrfitsReader();
	void check();
	void read_header();
	size_t read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading = false);
	MJD get_start_mjd_curfile(){return psf[idmap[ifile_cur]].primary.start_mjd;}
	double get_tsamp_curfile(){return psf[idmap[ifile_cur]].subint.tbin;}
	size_t get_count_curfile(){return ns_psfn;}
	size_t get_count(){return count;}
	size_t get_ifile(){return ifile_cur;}
	size_t get_ifile_ordered(){return idmap[ifile_cur];}
	void get_filterbank_template(Filterbank &filtem);

private:
	Integration it;
	Integration it8;
	size_t ntot;
	size_t count;
	size_t ns_psfn;
	std::vector<Psrfits> psf;
	long int ifile_cur;
	long int isubint_cur;
	long int isample_cur;
	bool update_file;
	bool update_subint;
};


#endif /* PSRFITSREADER_H */
