/*
 * integration.h
 *
 *  Created on: Feb 22, 2020
 *      Author: ypmen
 */

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

class Integration
{
public:
	enum Mode{SEARCH, FOLD};
	enum DataType{USHORT, SHORT, UINT1, UINT2, UINT4, UINT8, FLOAT};
public:
	Integration();
	Integration(const Integration &it);
	Integration & operator=(const Integration &it);
	~Integration();
	void free();
	void load_data(void *dat, int np, int nc, int nb);
	void load_frequencies(double *freq, int nc);
	void load_weights(float *wts, int nc);
	void resize(int np, int nc, int nb);
	bool to_char(Integration &it);
public:
	enum Mode mode;
	enum DataType dtype;

	double indexval;
	double folding_period;
	double tsubint;
	double offs_sub;
	double lst_sub;
	double ra_sub;
	double dec_sub;
	double glon_sub;
	double glat_sub;
	float fd_ang;
	float pos_ang;
	float par_ang;
	float tel_az;
	float tel_zen;
	float aux_dm;
	float aux_rm;

	int npol;
	int nchan;
	int nbin;
	int nsblk;
	int nbits;

	double * frequencies;
	float * weights;
	float * offsets;
	float * scales;
	void * data;
};



#endif /* INTEGRATION_H_ */
