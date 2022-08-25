/*
 * hdu.h
 *
 *  Created on: Feb 25, 2020
 *      Author: ypmen
 */

#ifndef HDU_H_
#define HDU_H_

#include <vector>
#include <string.h>
#include <fitsio.h>

#include "integration.h"
#include "mjd.h"

using namespace std;

class PrimaryHDU
{
public:
	PrimaryHDU();
	~PrimaryHDU();
	bool load(fitsfile *fptr);
	bool unload(fitsfile *fptr);
public:
	MJD start_mjd;
	char observer[FLEN_VALUE];
	char projid[FLEN_VALUE];
	char telesop[FLEN_VALUE];
	char ibeam[FLEN_VALUE];
	char obs_mode[FLEN_VALUE];
	char date_obs[FLEN_VALUE];
	double chan_dm;
	char src_name[FLEN_VALUE];
	char ra[FLEN_VALUE];
	char dec[FLEN_VALUE];
	char stt_crd1[FLEN_VALUE];
	char stt_crd2[FLEN_VALUE];
	char trk_mode[FLEN_VALUE];
	char stp_crd1[FLEN_VALUE];
	char stp_crd2[FLEN_VALUE];
	char fd_mode[FLEN_VALUE];
};

class PsrparamHDU
{
public:
	PsrparamHDU();
	PsrparamHDU(const string file);
	~PsrparamHDU();
	void load(const string file);
	bool load(fitsfile *fptr);
	bool unload(fitsfile *fptr);
public:
	vector<string> text;
};

class T2predictHDU
{
public:
	T2predictHDU();
	T2predictHDU(const string file);
	~T2predictHDU();
	void load(const string file);
	bool load(fitsfile *fptr);
	bool unload(fitsfile *fptr);
public:
	vector<string> text;
};

class SubintHDU
{
public:
	SubintHDU();
	~SubintHDU();
	void free();
	void resize(int ns);
	void load_data(float *profiles, int ns, int np, int nc, int nb);
	void load_data(unsigned char *profiles, int ns, int np, int nc, int nb);
	void load_frequencies(double *freq, int nc);
	void load_weights(float *wts, int nc);

	bool load(fitsfile *fptr);
	bool load_header(fitsfile *fptr);
	bool load_data(fitsfile *fptr);
	bool load_integration(fitsfile *fptr, int k);
	bool load_integration(fitsfile *fptr, int k, Integration &it);
	bool load_integration_data(fitsfile *fptr, int k, Integration &it);

	bool unload(fitsfile *fptr);
	bool unload_header(fitsfile *fptr);
	bool unload_data(fitsfile *fptr);
	bool unload_integration(fitsfile *fptr, int k);
	bool unload_integration(fitsfile *fptr, Integration &it);
public:
	int nbits;
	double tbin;
	int nsubint;
	int npol;
	int nchan;
	int nbin;

	int nsblk;
	long int nsuboffs;
	long int nstot;

	int nsamples;

	double zero_off;
	double dm;
	double rm;

	enum Integration::Mode mode;
	enum Integration::DataType dtype;

	Integration *integrations;
};

#endif /* HDU_H_ */
