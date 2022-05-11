/*
 * psrfits.cpp
 *
 *  Created on: Feb 20, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <type_traits>
#include <fitsio.h>

#include "psrfits.h"

using namespace std;

Psrfits::Psrfits()
{
	mode = Integration::FOLD;

	subint.mode = mode;
	fptr = NULL;
}

Psrfits::Psrfits(const string fname, int iomode)
{
	filename = fname;

	mode = Integration::FOLD;

	subint.mode = mode;
	fptr = NULL;
}

Psrfits::~Psrfits()
{
	if (fptr != NULL)
	{
		int status = 0;
		fits_close_file(fptr, &status);
		if (status != 0)
			fits_report_error(stderr, status);
		fptr = NULL;
	}
}

void Psrfits::close()
{
	subint.free();

	if (fptr != NULL)
	{
		int status = 0;
		fits_close_file(fptr, &status);
		if (status != 0)
			fits_report_error(stderr, status);
		fptr = NULL;
	}
}

void Psrfits::load_mode()
{
	if (strcmp(primary.obs_mode, "PSR") == 0)
		mode = Integration::FOLD;
	else
		mode = Integration::SEARCH;

	subint.mode = mode;
}

bool Psrfits::parse_template(const string temfile)
{
	if (!check_template(temfile))
	{
		cerr<<"Error: check template"<<endl;
		return false;
	}

	int status = 0;
	fits_create_file (&fptr, (filename).c_str(), &status);
	if (status)
	{
		cerr<<"Error: can not open file "<<filename<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_execute_template(fptr, const_cast<char*>(temfile.c_str()), &status);
	if (status)
	{
		cerr<<"Error: can not execute template file "<<temfile<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

//from psrchive
bool Psrfits::check_template(const string temfile)
{
	FILE* fptr = fopen(temfile.c_str(), "r");
	if (!fptr)
	{
		cerr<<"Error: can not open file "<<temfile<<endl;
		return false;
	}

	char templt[FLEN_CARD*2];
	char card[FLEN_CARD*2];

	int status = 0;

	while (fgets (templt, FLEN_CARD*2, fptr))
	{
		// CFITSIO User's Reference Guide
		// 11.1 Detailed Template Line Format

		/* "Any template line that begins with the pound '#' character is
		   ignored by the template parser and may be use to insert
		   comments into the template file itself." */

		if (templt[0] == '#')
		  continue;

		char* newline = strchr (templt, '\n');
		if (newline)
		  *newline = '\0';

		int keytype = 0;
		fits_parse_template (templt, card, &keytype, &status);
		if (status)
		{
			cerr<<"Error: template '"<<templt<<"'"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (strlen(templt) >= FLEN_CARD)
		{
			cerr<<"Error: card '"<<templt<<"' out of length"<<endl;
			return false;
		}

		fits_test_record (card, &status);
		if (status)
		{
			cerr<<"Error: card '"<<card<<"'"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		char keyname[FLEN_CARD];
		int keylength = 0;
		fits_get_keyname (card, keyname, &keylength, &status);
		if (status)
		{
			cerr<<"Error: card '"<<card<<"'"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (keyname[strlen(keyname)-1]=='#')
		  keyname[strlen(keyname)-1] = '\0'; // ignore auto-indexing marks

		fits_test_keyword (keyname, &status);
		if (status)
		{
			cerr<<"Error: keyname '"<<keyname<<"'"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	  }

	  fclose (fptr);

	  return true;
}

bool Psrfits::open(int iomode)
{
	int status = 0;

	fits_open_file(&fptr, filename.c_str(), iomode, &status);
	if (status)
	{
		cerr<<"Error: can not open file '"<<filename<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool Psrfits::remove_hdu(const string hduname)
{
	int status=0;

	char extname[FLEN_VALUE];
	strcpy(extname, hduname.c_str());
	fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to hdu "<<hduname<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_delete_hdu(fptr, NULL, &status);
	if(status)
	{
		cerr<<"Error: can not remove hdu "<<hduname<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}
