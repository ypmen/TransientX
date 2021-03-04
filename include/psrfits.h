/*
 * psrfits.h
 *
 *  Created on: Feb 19, 2020
 *      Author: ypmen
 */

#ifndef PSRFITS_H_
#define PSRFITS_H_

#include <vector>
#include <string.h>
#include <fitsio.h>

#include "hdu.h"

using namespace std;

class Psrfits
{
public:
	Psrfits();
	Psrfits(const string fname, int iomode=0);
	~Psrfits();
	bool open(int iomode=0);
	void close();
	void load_mode();
	static bool check_template(const string temfile);
	bool parse_template(const string temfile);
	bool remove_hdu(const string hduname);
public:
	string filename;

	enum Integration::Mode mode;

	PrimaryHDU primary;
	PsrparamHDU psrparam;
	T2predictHDU t2predict;
	SubintHDU subint;
	fitsfile *fptr;
};


#endif /* PSRFITS_H_ */
