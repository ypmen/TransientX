/*
 * hdu.cpp
 *
 *  Created on: Feb 25, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <fstream>
#include <fitsio.h>

#include "psrfits.h"
#include "hdu.h"

using namespace std;

PrimaryHDU::PrimaryHDU()
{
	strcpy(observer, "");
	strcpy(projid, "");
	strcpy(telesop, "");
	strcpy(ibeam, "");
	strcpy(obs_mode, "");
	strcpy(date_obs, "");
	strcpy(src_name, "");
	strcpy(ra, "");
	strcpy(dec, "");
	strcpy(stt_crd1, "");
	strcpy(stt_crd2, "");
	strcpy(trk_mode, "");
	strcpy(stp_crd1, "");
	strcpy(stp_crd2, "");
	strcpy(fd_mode, "");
	chan_dm = 0.;
}

PrimaryHDU::~PrimaryHDU()
{
}

bool PrimaryHDU::load(fitsfile *fptr)
{
	bool res = true;

	int status = 0;

	fits_movabs_hdu(fptr, 1, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not move to PRIMARY"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_key(fptr, TSTRING, "OBSERVER", observer, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read OBSERVER"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "PROJID", projid, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read PROJID"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "TELESCOP", telesop, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read TELESCOP"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "IBEAM", ibeam, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read IBEAM"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "OBS_MODE", obs_mode, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "DATE-OBS", date_obs, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read DATE-OBS"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TDOUBLE, "CHAN_DM", &chan_dm, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read CHAN_DM"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "SRC_NAME", src_name, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read SRC_NAME"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "RA", ra, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read RA"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "DEC", dec, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read DEC"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "STT_CRD1", stt_crd1, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "STT_CRD2", stt_crd2, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "TRK_MODE", trk_mode, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read TRK_MODE"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "STP_CRD1", stp_crd1, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read STP_CRD1"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "STP_CRD2", stp_crd2, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read STP_CRD2"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TSTRING, "FD_MODE", fd_mode, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read FD_MODE"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TLONG, "STT_IMJD", &(start_mjd.stt_imjd), NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read STT_IMJD"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TLONG, "STT_SMJD", &(start_mjd.stt_smjd), NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read STT_SMJD"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	fits_read_key(fptr, TDOUBLE, "STT_OFFS", &(start_mjd.stt_offs), NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read STT_OFFS"<<endl;
		fits_report_error(stderr, status);
		status = 0;
		res = false;
	}

	return res;
}

bool PrimaryHDU::unload(fitsfile *fptr)
{
	int status = 0;

	fits_movabs_hdu(fptr, 1, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not move to PRIMARY"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_write_date(fptr, &status);
	if (status)
	{
		cerr<<"Error: can not write date"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "OBSERVER", observer, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write OBSERVER"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "PROJID", projid, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write PROJID"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "TELESCOP", telesop, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write TELESCOP"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "IBEAM", ibeam, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write IBEAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "OBS_MODE", obs_mode, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "DATE-OBS", date_obs, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write DATE-OBS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TDOUBLE, "CHAN_DM", &chan_dm, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write CHAN_DM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "SRC_NAME", src_name, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write SRC_NAME"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "RA", ra, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write RA"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "DEC", dec, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write DEC"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "STT_CRD1", stt_crd1, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "STT_CRD2", stt_crd2, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write OBS_MODE"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "TRK_MODE", trk_mode, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write TRK_MODE"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "STP_CRD1", stp_crd1, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write STP_CRD1"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "STP_CRD2", stp_crd2, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write STP_CRD2"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TSTRING, "FD_MODE", fd_mode, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write FD_MODE"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	start_mjd.format();

	fits_update_key(fptr, TLONG, "STT_IMJD", &(start_mjd.stt_imjd), NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set STT_IMJD"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TLONG, "STT_SMJD", &(start_mjd.stt_smjd), NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set STT_SMJD"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TDOUBLE, "STT_OFFS", &(start_mjd.stt_offs), NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set STT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

PsrparamHDU::PsrparamHDU()
{
	text.reserve(128);
}
PsrparamHDU::PsrparamHDU(const string file)
{
	text.reserve(128);
	load(file);
}
PsrparamHDU::~PsrparamHDU()
{
}

void PsrparamHDU::load(const string file)
{
	ifstream infile(file);
	string line;
	while (getline(infile, line))
	{
		text.push_back(line);
	}
	infile.close();
}

bool PsrparamHDU::load(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "PSRPARAM", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to PSRPARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	long int nrows = 0;
	fits_get_num_rows(fptr, &nrows, &status);
	if (status)
	{
		cerr<<"Error: can not read nrows"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (nrows == 0) return false;

	int colnum = 0;

	fits_get_colnum (fptr, CASEINSEN, "PARAM", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of PARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	int typecode;
	long int repeat = 0;
	long int width = 0;
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	if (status)
	{
		cerr<<"Error: can not read coltype of PARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	char *line = new char [repeat+1];
	for (long int l=0; l<nrows; l++)
	{
		fits_read_col(fptr, TSTRING, colnum, l+1, 1, 1, 0, &line, 0, &status);
		text.push_back(line);
	}
	delete [] line;

	if (status)
	{
		cerr<<"Error: can not read column PARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool PsrparamHDU::unload(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "PSRPARAM", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to PSRPARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	int colnum = 0;

	fits_get_colnum (fptr, CASEINSEN, "PARAM", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of PARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<text.size(); l++)
	{
		fits_write_col(fptr, TSTRING, colnum, l+1, 1, 1, &(text[l]), &status);
	}

	if (status)
	{
		cerr<<"Error: can not set column PARAM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

T2predictHDU::T2predictHDU()
{
	text.reserve(128);
}

T2predictHDU::T2predictHDU(const string file)
{
	text.reserve(128);
	load(file);
}

T2predictHDU::~T2predictHDU()
{
}

void T2predictHDU::load(const string file)
{
	ifstream infile(file);
	string line;
	while (getline(infile, line))
	{
		text.push_back(line);
	}
	infile.close();
}

bool T2predictHDU::load(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "T2PREDICT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to T2PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	long int nrows = 0;
	fits_get_num_rows(fptr, &nrows, &status);
	if (status)
	{
		cerr<<"Error: can not read nrows"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (nrows == 0) return false;

	int colnum = 0;

	fits_get_colnum (fptr, CASEINSEN, "PREDICT", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	int typecode;
	long int repeat = 0;
	long int width = 0;
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	if (status)
	{
		cerr<<"Error: can not read coltype of PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	char *line = new char [repeat+1];
	for (long int l=0; l<nrows; l++)
	{
		fits_read_col(fptr, TSTRING, colnum, l+1, 1, 1, 0, &line, 0, &status);
		text.push_back(line);
	}
	delete [] line;

	if (status)
	{
		cerr<<"Error: can not set column PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool T2predictHDU::unload(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "T2PREDICT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to T2PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	int colnum = 0;

	fits_get_colnum (fptr, CASEINSEN, "PREDICT", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<text.size(); l++)
	{
		fits_write_col(fptr, TSTRING, colnum, l+1, 1, 1, &(text[l]), &status);
	}

	if (status)
	{
		cerr<<"Error: can not set column PREDICT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

SubintHDU::SubintHDU()
{
	nbits = 1;

	tbin = 0.;
	nsubint = 0;
	npol = 0;
	nchan = 0;
	nbin = 0;

	nsblk = 1;
	nsuboffs = 0;
	nstot = 1;

	nsamples = 0;

	zero_off = 0.;
	dm = 0.;
	rm = 0.;

	mode = Integration::FOLD;
	dtype = Integration::SHORT;

	integrations = NULL;
}

SubintHDU::~SubintHDU()
{
	if (integrations != NULL)
	{
		delete [] integrations;
		integrations = NULL;
	}
}

void SubintHDU::free()
{
	if (integrations != NULL)
	{
		delete [] integrations;
		integrations = NULL;
	}
}

void SubintHDU::resize(int ns)
{
	nsubint = ns;
	if (mode == Integration::FOLD)
	{
		if (integrations != NULL)
			delete [] integrations;

		integrations = new Integration [nsubint];
		for (long int l=0; l<nsubint; l++)
		{
			integrations[l].mode = mode;
			integrations[l].dtype = dtype;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		if (integrations != NULL) delete [] integrations;
		integrations = new Integration [nsubint];
		for (long int l=0; l<nsubint; l++)
		{
			integrations[l].mode = mode;
			integrations[l].dtype = dtype;
		}
	}
}

void SubintHDU::load_data(float *profiles, int ns, int np, int nc, int nb)
{
	resize(ns);
	npol = np;
	nchan = nc;
	nbin = nb;

	float *pro = profiles;
	for (long int l=0; l<nsubint; l++)
	{
		integrations[l].load_data(pro, npol, nchan, nbin);
		pro += npol*nchan*nbin;
	}
	pro = NULL;
}

void SubintHDU::load_data(unsigned char *profiles, int ns, int np, int nc, int nb)
{
	resize(ns);
	npol = np;
	nchan = nc;
	nsblk = nb;

	unsigned char *pro = profiles;
	for (long int l=0; l<nsubint; l++)
	{
		integrations[l].load_data(pro, npol, nchan, nsblk);
		pro += nsblk*npol*nchan*nbits/8;
	}
	pro = NULL;

}

void SubintHDU::load_frequencies(double *freq, int nc)
{
	int nch = nc<nchan ? nc:nchan;

	for (long int l=0; l<nsubint; l++)
	{
		integrations[l].load_frequencies(freq, nch);
	}
}

void SubintHDU::load_weights(float *wts, int nc)
{
	int nch = nc<nchan ? nc:nchan;

	for (long int l=0; l<nsubint; l++)
	{
		integrations[l].load_weights(wts, nch);
	}
}

bool SubintHDU::load(fitsfile *fptr)
{
	if (load_header(fptr) and load_data(fptr))
		return true;

	return false;
}

bool SubintHDU::load_header(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	long int nrows = 0;
	fits_get_num_rows(fptr, &nrows, &status);
	if (status)
	{
		cerr<<"Error: can not read nrows"<<endl;
		fits_report_error(stderr, status);
		return false;
	}
	nsubint = nrows;

	fits_read_key(fptr, TINT, "NPOL", &npol, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not read NPOL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_key(fptr, TINT, "NCHAN", &nchan, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not read NCHAN"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_key(fptr, TINT, "NBIN", &nbin, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not read NBIN"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_key(fptr, TINT, "NBITS", &nbits, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not read nbits"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (mode == Integration::SEARCH)
	{
		switch (nbits)
		{
		case 1: dtype=Integration::UINT1; break;
		case 2: dtype=Integration::UINT2; break;
		case 4: dtype=Integration::UINT4; break;
		case 8: dtype=Integration::UINT8; break;
		case 32: dtype=Integration::FLOAT; break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
	}

	fits_read_key(fptr, TLONG, "NSUBOFFS", &nsuboffs, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read NSUBOFFS"<<endl;
		status = 0;
	}

	fits_read_key(fptr, TINT, "NSBLK", &nsblk, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not read NSBLK"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_key(fptr, TLONG, "NSTOT", &nstot, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read NSTOT"<<endl;
		status = 0;
	}

	fits_read_key(fptr, TDOUBLE, "TBIN", &tbin, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read TBIN"<<endl;
		status = 0;
	}

	fits_read_key(fptr, TDOUBLE, "ZERO_OFF", &zero_off, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read ZERO_OFF"<<endl;
		status = 0;
	}

	fits_read_key(fptr, TDOUBLE, "DM", &dm, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read DM"<<endl;
		status = 0;
	}

	fits_read_key(fptr, TDOUBLE, "RM", &rm, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not read RM"<<endl;
		status = 0;
	}

	if (nstot > 1)
		nsamples = nstot;
	else
		nsamples = nsubint*nsblk;

	resize(nsubint);

	return true;
}

bool SubintHDU::load_data(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status == 0)
		{
			for (long int l=0; l<nsubint; l++)
			{
				fits_read_col(fptr, TDOUBLE, colnum, l+1, 1, 1, 0, &(integrations[l].folding_period), 0, &status);
			}
			if (status)
			{
				cerr<<"Error: can not read PERIOD"<<endl;
				fits_report_error(stderr, status);
				return false;
			}
		}
		else
		{
			cerr<<"Warning: can not read colnum of PERIOD"<<endl;
			status = 0;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status == 0)
	{
		for (long int l=0; l<nsubint; l++)
		{
			fits_read_col(fptr, TDOUBLE, colnum, l+1, 1, 1, 0, &(integrations[l].tsubint), 0, &status);
		}
		if (status)
		{
			cerr<<"Error: can not read TSUBINT"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of TSUBINT"<<endl;
		status = 0;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status == 0)
	{
		for (long int l=0; l<nsubint; l++)
		{
			fits_read_col(fptr, TDOUBLE, colnum, l+1, 1, 1, 0, &(integrations[l].offs_sub), 0, &status);
		}
		if (status)
		{
			cerr<<"Error: can not read OFFS_SUB"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of OFFS_SUB"<<endl;
		status = 0;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		for (long int l=0; l<nsubint; l++)
		{
			integrations[l].mode = mode;
			integrations[l].dtype = dtype;
			integrations[l].resize(npol, nchan, nbin);

			fits_read_col(fptr, TSHORT, colnum, l+1, 1, npol*nchan*nbin, 0, integrations[l].data, 0, &status);
		}
		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{

		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		for (long int l=0; l<nsubint; l++)
		{
			integrations[l].mode = mode;
			integrations[l].resize(npol, nchan, nsblk);

			switch (dtype)
			{
			case Integration::UINT1: fits_read_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[l].data, 0, &status); break;
			case Integration::UINT2: fits_read_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[l].data, 0, &status); break;
			case Integration::UINT4: fits_read_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[l].data, 0, &status); break;
			case Integration::UINT8: fits_read_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[l].data, 0, &status); break;
			case Integration::FLOAT: fits_read_col(fptr, TFLOAT, colnum, l+1, 1, nsblk*npol*nchan, 0, integrations[l].data, 0, &status); break;
			default: cerr<<"Error: data type not supported"<<endl; break;
			}
		}
		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_read_col(fptr, TDOUBLE, colnum, l+1, 1, nchan, 0, integrations[l].frequencies, 0, &status);
	}
	if (status)
	{
		cerr<<"Error: can not read DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_read_col(fptr, TFLOAT, colnum, l+1, 1, nchan, 0, integrations[l].weights, 0, &status);
	}
	if (status)
	{
		cerr<<"Error: can not read DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_read_col(fptr, TFLOAT, colnum, l+1, 1, npol*nchan, 0, integrations[l].scales, 0, &status);
	}
	if (status)
	{
		cerr<<"Error: can not read DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_read_col(fptr, TFLOAT, colnum, l+1, 1, npol*nchan, 0, integrations[l].offsets, 0, &status);
	}
	if (status)
	{
		cerr<<"Error: can not read DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::load_integration(fitsfile *fptr, int k)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status == 0)
		{
			fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(integrations[k].folding_period), 0, &status);
			if (status)
			{
				cerr<<"Error: can not read PERIOD"<<endl;
				fits_report_error(stderr, status);
				return false;
			}
		}
		else
		{
			cerr<<"Warning: can not read colnum of PERIOD"<<endl;
			status = 0;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status == 0)
	{
		fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(integrations[k].tsubint), 0, &status);
		if (status)
		{
			cerr<<"Error: can not read TSUBINT"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of TSUBINT"<<endl;
		status = 0;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status == 0)
	{
		fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(integrations[k].offs_sub), 0, &status);
		if (status)
		{
			cerr<<"Error: can not read OFFS_SUB"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of OFFS_SUB"<<endl;
		status = 0;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		integrations[k].mode = mode;
		integrations[k].dtype = dtype;
		integrations[k].resize(npol, nchan, nbin);

		fits_read_col(fptr, TSHORT, colnum, k+1, 1, npol*nchan*nbin, 0, integrations[k].data, 0, &status);

		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		integrations[k].mode = mode;
		integrations[k].dtype = dtype;
		integrations[k].resize(npol, nchan, nsblk);

		switch (dtype)
		{
		case Integration::UINT1: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[k].data, 0, &status); break;
		case Integration::UINT2: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[k].data, 0, &status); break;
		case Integration::UINT4: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[k].data, 0, &status); break;
		case Integration::UINT8: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, integrations[k].data, 0, &status); break;
		case Integration::FLOAT: fits_read_col(fptr, TFLOAT, colnum, k+1, 1, nsblk*npol*nchan, 0, integrations[k].data, 0, &status); break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, nchan, 0, integrations[k].frequencies, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, nchan, 0, integrations[k].weights, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, 0, integrations[k].scales, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, 0, integrations[k].offsets, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::load_integration(fitsfile *fptr, int k, Integration &it)
{

	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status == 0)
		{
			fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(it.folding_period), 0, &status);
			if (status)
			{
				cerr<<"Error: can not read PERIOD"<<endl;
				fits_report_error(stderr, status);
				return false;
			}
		}
		else
		{
			cerr<<"Warning: can not read colnum of PERIOD"<<endl;
			status = 0;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status == 0)
	{
		fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(it.tsubint), 0, &status);
		if (status)
		{
			cerr<<"Error: can not read TSUBINT"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of TSUBINT"<<endl;
		status = 0;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status == 0)
	{
		fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, 1, 0, &(it.offs_sub), 0, &status);
		if (status)
		{
			cerr<<"Error: can not read OFFS_SUB"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else
	{
		cerr<<"Warning: can not read colnum of OFFS_SUB"<<endl;
		status = 0;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		it.mode = mode;
		it.dtype = dtype;
		it.resize(npol, nchan, nbin);

		fits_read_col(fptr, TSHORT, colnum, k+1, 1, npol*nchan*nbin, 0, it.data, 0, &status);

		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		it.mode = mode;
		it.dtype = dtype;
		it.resize(npol, nchan, nsblk);

		switch (dtype)
		{
		case Integration::UINT1: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT2: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT4: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT8: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::FLOAT: fits_read_col(fptr, TFLOAT, colnum, k+1, 1, nsblk*npol*nchan, 0, it.data, 0, &status); break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TDOUBLE, colnum, k+1, 1, nchan, 0, it.frequencies, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, nchan, 0, it.weights, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, 0, it.scales, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_read_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, 0, it.offsets, 0, &status);
	if (status)
	{
		cerr<<"Error: can not read DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::load_integration_data(fitsfile *fptr, int k, Integration &it)
{

	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		it.mode = mode;
		it.dtype = dtype;
		it.resize(npol, nchan, nbin);

		fits_read_col(fptr, TSHORT, colnum, k+1, 1, npol*nchan*nbin, 0, it.data, 0, &status);

		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		it.mode = mode;
		it.dtype = dtype;
		it.resize(npol, nchan, nsblk);

		switch (dtype)
		{
		case Integration::UINT1: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT2: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT4: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::UINT8: fits_read_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, 0, it.data, 0, &status); break;
		case Integration::FLOAT: fits_read_col(fptr, TFLOAT, colnum, k+1, 1, nsblk*npol*nchan, 0, it.data, 0, &status); break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
		if (status)
		{
			cerr<<"Error: can not read DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	return true;
}

bool SubintHDU::unload(fitsfile *fptr)
{
	if (unload_header(fptr) and unload_data(fptr))
		return true;
	return false;
}

bool SubintHDU::unload_header(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TINT, "NPOL", &npol, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set NPOL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TDOUBLE, "TBIN", &tbin, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set TBIN"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TINT, "NBIN", &nbin, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set NBIN"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	switch (dtype)
	{
	case Integration::UINT1: nbits=1; break;
	case Integration::UINT2: nbits=2; break;
	case Integration::UINT4: nbits=4; break;
	case Integration::UINT8: nbits=8; break;
	case Integration::FLOAT: nbits=32; break;
	case Integration::SHORT: nbits=1; break;
	case Integration::USHORT: nbits=1; break;
	default: cerr<<"Error: data type not supported"<<endl; break;
	}
	fits_update_key(fptr, TINT, "NBITS", &nbits, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write nbits"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TINT, "NCHAN", &nchan, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set NCHAN"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TLONG, "NSUBOFFS", &nsuboffs, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not write NSUBOFFS"<<endl;
		status = 0;
	}

	fits_update_key(fptr, TINT, "NSBLK", &nsblk, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not write NSBLK"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TLONG, "NSTOT", &nstot, NULL, &status);
	if (status)
	{
		cerr<<"Warning: can not write NSTOT"<<endl;
		status = 0;
	}

	fits_update_key(fptr, TDOUBLE, "DM", &dm, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set DM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_update_key(fptr, TDOUBLE, "RM", &rm, NULL, &status);
	if (status)
	{
		cerr<<"Error: can not set RM"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::unload_data(fitsfile *fptr)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		for (long int l=0; l<nsubint; l++)
		{
			fits_write_col(fptr, TDOUBLE, colnum, l+1, 1, 1, &(integrations[l].folding_period), &status);
		}
		if (status)
		{
			cerr<<"Error: can not set PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TDOUBLE, colnum, l+1, 1, 1, &(integrations[l].tsubint), &status);
	}
	if (status)
	{
		cerr<<"Error: can not set TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TDOUBLE, colnum, l+1, 1, 1, &(integrations[l].offs_sub), &status);
	}
	if (status)
	{
		cerr<<"Error: can not set OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		fits_modify_vector_len(fptr, colnum, npol*nchan*nbin, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		int naxis = 3;
		long int naxes[3] = {nbin, nchan, npol};
		fits_write_tdim(fptr, colnum, naxis, naxes, &status);

		for (long int l=0; l<nsubint; l++)
		{
			fits_write_col(fptr, TSHORT, colnum, l+1, 1, npol*nchan*nbin, integrations[l].data, &status);
		}
		if (status)
		{
			cerr<<"Error: can not set DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (nbits != 32)
			fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan*nbits/8, &status);
		else
			fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan, &status);

		if (status)
		{
			cerr<<"Error: can not resize column DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		int naxis = 3;
		long int naxes[3] = {nchan, npol, nsblk};
		fits_write_tdim(fptr, colnum, naxis, naxes, &status);

		for (long int l=0; l<nsubint; l++)
		{
			switch (dtype)
			{
			case Integration::UINT1: fits_write_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, integrations[l].data, &status); break;
			case Integration::UINT2: fits_write_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, integrations[l].data, &status); break;
			case Integration::UINT4: fits_write_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, integrations[l].data, &status); break;
			case Integration::UINT8: fits_write_col(fptr, TBYTE, colnum, l+1, 1, nsblk*npol*nchan*nbits/8, integrations[l].data, &status); break;
			case Integration::FLOAT: fits_write_col(fptr, TFLOAT, colnum, l+1, 1, nsblk*npol*nchan, integrations[l].data, &status); break;
			default: cerr<<"Error: data type not supported"<<endl; break;
			}
		}
		if (status)
		{
			cerr<<"Error: can not write DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_modify_vector_len(fptr, colnum, nchan, &status);
	if (status)
	{
		cerr<<"Error: can not resize column DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TDOUBLE, colnum, l+1, 1, nchan, integrations[l].frequencies, &status);
	}
	if (status)
	{
		cerr<<"Error: can not set DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_modify_vector_len(fptr, colnum, nchan, &status);
	if (status)
	{
		cerr<<"Error: can not resize column DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TFLOAT, colnum, l+1, 1, nchan, integrations[l].weights, &status);
	}
	if (status)
	{
		cerr<<"Error: can not set DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
	if (status)
	{
		cerr<<"Error: can not resize column DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TFLOAT, colnum, l+1, 1, npol*nchan, integrations[l].scales, &status);
	}
	if (status)
	{
		cerr<<"Error: can not set DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
	if (status)
	{
		cerr<<"Error: can not resize column DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	for (long int l=0; l<nsubint; l++)
	{
		fits_write_col(fptr, TFLOAT, colnum, l+1, 1, npol*nchan, integrations[l].offsets, &status);
	}
	if (status)
	{
		cerr<<"Error: can not set DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::unload_integration(fitsfile *fptr, int k)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(integrations[k].folding_period), &status);
		if (status)
		{
			cerr<<"Error: can not set PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(integrations[k].tsubint), &status);
	if (status)
	{
		cerr<<"Error: can not set TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(integrations[k].offs_sub), &status);
	if (status)
	{
		cerr<<"Error: can not set OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (k==0)
		{
			fits_modify_vector_len(fptr, colnum, npol*nchan*nbin, &status);
			if (status)
			{
				cerr<<"Error: can not resize column DATA"<<endl;
				fits_report_error(stderr, status);
				return false;
			}

			int naxis = 3;
			long int naxes[3] = {nbin, nchan, npol};
			fits_write_tdim(fptr, colnum, naxis, naxes, &status);
		}

		fits_write_col(fptr, TSHORT, colnum, k+1, 1, npol*nchan*nbin, integrations[k].data, &status);
		if (status)
		{
			cerr<<"Error: can not set DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (k==0)
		{
			if (nbits != 32)
				fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan*nbits/8, &status);
			else
				fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan, &status);

			if (status)
			{
				cerr<<"Error: can not resize column DATA"<<endl;
				fits_report_error(stderr, status);
				return false;
			}

			int naxis = 3;
			long int naxes[3] = {nchan, npol, nsblk};
			fits_write_tdim(fptr, colnum, naxis, naxes, &status);
		}

		switch (dtype)
		{
		case Integration::UINT1: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, integrations[k].data, &status); break;
		case Integration::UINT2: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, integrations[k].data, &status); break;
		case Integration::UINT4: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, integrations[k].data, &status); break;
		case Integration::UINT8: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, integrations[k].data, &status); break;
		case Integration::FLOAT: fits_write_col(fptr, TFLOAT, colnum, k+1, 1, nsblk*npol*nchan, integrations[k].data, &status); break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
		if (status)
		{
			cerr<<"Error: can not write DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_FREQ"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, nchan, integrations[k].frequencies, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_WTS"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, nchan, integrations[k].weights, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_SCL"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, integrations[k].scales, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_OFFS"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, integrations[k].offsets, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}

bool SubintHDU::unload_integration(fitsfile *fptr, Integration &it)
{
	int status = 0;

	fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status);
	if (status)
	{
		cerr<<"Error: can not move to SUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	long int k = 0;
	fits_get_num_rows(fptr, &k, &status);
	if (status)
	{
		cerr<<"Error: can not read nrows"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//table
	int colnum = 0;

	if (mode == Integration::FOLD)
	{
		//PERIOD
		fits_get_colnum(fptr, CASEINSEN, "PERIOD", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(it.folding_period), &status);
		if (status)
		{
			cerr<<"Error: can not set PERIOD"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//TSUBINT
	fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(it.tsubint), &status);
	if (status)
	{
		cerr<<"Error: can not set TSUBINT"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//OFFS_SUB
	fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, 1, &(it.offs_sub), &status);
	if (status)
	{
		cerr<<"Error: can not set OFFS_SUB"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DATA
	if (mode == Integration::FOLD)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (k==0)
		{
			fits_modify_vector_len(fptr, colnum, npol*nchan*nbin, &status);
			if (status)
			{
				cerr<<"Error: can not resize column DATA"<<endl;
				fits_report_error(stderr, status);
				return false;
			}

			int naxis = 3;
			long int naxes[3] = {nbin, nchan, npol};
			fits_write_tdim(fptr, colnum, naxis, naxes, &status);
		}

		fits_write_col(fptr, TSHORT, colnum, k+1, 1, npol*nchan*nbin, it.data, &status);
		if (status)
		{
			cerr<<"Error: can not set DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}
	else if (mode == Integration::SEARCH)
	{
		fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		if (status)
		{
			cerr<<"Error: can not read column number of DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}

		if (k==0)
		{
			if (nbits != 32)
				fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan*nbits/8, &status);
			else
				fits_modify_vector_len(fptr, colnum, nsblk*npol*nchan, &status);
			
			if (status)
			{
				cerr<<"Error: can not resize column DATA"<<endl;
				fits_report_error(stderr, status);
				return false;
			}

			int naxis = 3;
			long int naxes[3] = {nchan, npol, nsblk};
			if (nbits != 32)
				naxes[2] /= 8/nbits;
			fits_write_tdim(fptr, colnum, naxis, naxes, &status);
		}

		switch (dtype)
		{
		case Integration::UINT1: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, it.data, &status); break;
		case Integration::UINT2: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, it.data, &status); break;
		case Integration::UINT4: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, it.data, &status); break;
		case Integration::UINT8: fits_write_col(fptr, TBYTE, colnum, k+1, 1, nsblk*npol*nchan*nbits/8, it.data, &status); break;
		case Integration::FLOAT: fits_write_col(fptr, TFLOAT, colnum, k+1, 1, nsblk*npol*nchan, it.data, &status); break;
		default: cerr<<"Error: data type not supported"<<endl; break;
		}
		if (status)
		{
			cerr<<"Error: can not write DATA"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	//DAT_FREQ
	fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_FREQ"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TDOUBLE, colnum, k+1, 1, nchan, it.frequencies, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_FREQ"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_WTS
	fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_WTS"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, nchan, it.weights, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_WTS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_SCL
	fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_SCL"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, it.scales, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_SCL"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	//DAT_OFFS
	fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &colnum, &status);
	if (status)
	{
		cerr<<"Error: can not read column number of DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	if (k==0)
	{
		fits_modify_vector_len(fptr, colnum, npol*nchan, &status);
		if (status)
		{
			cerr<<"Error: can not resize column DAT_OFFS"<<endl;
			fits_report_error(stderr, status);
			return false;
		}
	}

	fits_write_col(fptr, TFLOAT, colnum, k+1, 1, npol*nchan, it.offsets, &status);
	if (status)
	{
		cerr<<"Error: can not set DAT_OFFS"<<endl;
		fits_report_error(stderr, status);
		return false;
	}

	return true;
}
