/*
 * integration.cpp
 *
 *  Created on: Feb 22, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <cstring>
#include <assert.h>

#include "psrfits.h"
#include "hdu.h"

#include "integration.h"

using namespace std;

Integration::Integration()
{
	mode = FOLD;
	dtype = SHORT;

	indexval = 0.;
	folding_period = 0.;
	tsubint = 0.;
	offs_sub = 0.;
	lst_sub = 0.;
	ra_sub = 0.;
	dec_sub = 0.;
	glon_sub = 0.;
	glat_sub = 0.;
	fd_ang = 0.;
	pos_ang = 0.;
	par_ang = 0.;
	tel_az = 0.;
	tel_zen = 0.;
	aux_dm = 0.;
	aux_rm = 0.;

	npol = 0;
	nchan = 0;
	nbin = 0;
	nsblk = 1;
	nbits = 1;

	frequencies = NULL;
	weights = NULL;
	offsets = NULL;
	scales = NULL;
	data = NULL;
}

Integration::Integration(const Integration &it)
{
	mode = it.mode;
	dtype = it.dtype;

	indexval = it.indexval;
	folding_period = it.folding_period;
	tsubint = it.tsubint;
	offs_sub =it. offs_sub;
	lst_sub = it.lst_sub;
	ra_sub = it.ra_sub;
	dec_sub = it.dec_sub;
	glon_sub = it.glon_sub;
	glat_sub = it.glat_sub;
	fd_ang = it.fd_ang;
	pos_ang = it.pos_ang;
	par_ang = it.par_ang;
	tel_az = it.tel_az;
	tel_zen = it.tel_zen;
	aux_dm = it.aux_dm;
	aux_rm = it.aux_rm;

	npol = it.npol;
	nchan = it.nchan;
	nbin = it.nbin;
	nsblk = it.nsblk;
	nbits = it.nbits;

	if (it.frequencies != NULL)
	{
		if (frequencies != NULL) delete [] frequencies;
		frequencies = new double [nchan];
		std::memcpy(frequencies, it.frequencies, sizeof(double) * nchan);
	}
	else
	{
		frequencies = NULL;
	}

	if (it.weights != NULL)
	{
		if (weights != NULL) delete [] weights;
		weights = new float [nchan];
		std::memcpy(weights, it.weights, sizeof(float) * nchan);
	}
	else
	{
		weights = NULL;
	}

	if (it.offsets != NULL)
	{
		if (offsets != NULL) delete [] offsets;
		offsets = new float [npol * nchan];
		std::memcpy(offsets, it.offsets, sizeof(float) * npol * nchan);
	}
	else
	{
		offsets = NULL;
	}

	if (it.scales != NULL)
	{
		if (scales != NULL) delete [] scales;
		scales = new float [npol * nchan];
		std::memcpy(scales, it.scales, sizeof(float) * npol * nchan);
	}
	else
	{
		scales = NULL;
	}

	if (it.data != NULL)
	{
		if (data != NULL) delete [] data;
		if (mode == FOLD)
		{
			data = new short [npol * nchan * nbin];
			std::memcpy(data, it.data, sizeof(short) * npol * nchan * nbin);
		}
		else if (mode == SEARCH)
		{
			switch (dtype)
			{
			case UINT1: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT2: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT4: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT8: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case FLOAT: data = new float [nsblk * npol * nchan]; std::memcpy(data, it.data, sizeof(float) * nsblk * npol * nchan); break;
			default: cerr<<"Error: data type not support"<<endl; break;
			}
		}
	}
	else
	{
		data = NULL;
	}
}

Integration & Integration::operator=(const Integration &it)
{
	mode = it.mode;
	dtype = it.dtype;

	indexval = it.indexval;
	folding_period = it.folding_period;
	tsubint = it.tsubint;
	offs_sub =it. offs_sub;
	lst_sub = it.lst_sub;
	ra_sub = it.ra_sub;
	dec_sub = it.dec_sub;
	glon_sub = it.glon_sub;
	glat_sub = it.glat_sub;
	fd_ang = it.fd_ang;
	pos_ang = it.pos_ang;
	par_ang = it.par_ang;
	tel_az = it.tel_az;
	tel_zen = it.tel_zen;
	aux_dm = it.aux_dm;
	aux_rm = it.aux_rm;

	npol = it.npol;
	nchan = it.nchan;
	nbin = it.nbin;
	nsblk = it.nsblk;
	nbits = it.nbits;

	if (it.frequencies != NULL)
	{
		if (frequencies != NULL) delete [] frequencies;
		frequencies = new double [nchan];
		std::memcpy(frequencies, it.frequencies, sizeof(double) * nchan);
	}
	else
	{
		frequencies = NULL;
	}

	if (it.weights != NULL)
	{
		if (weights != NULL) delete [] weights;
		weights = new float [nchan];
		std::memcpy(weights, it.weights, sizeof(float) * nchan);
	}
	else
	{
		weights = NULL;
	}

	if (it.offsets != NULL)
	{
		if (offsets != NULL) delete [] offsets;
		offsets = new float [npol * nchan];
		std::memcpy(offsets, it.offsets, sizeof(float) * npol * nchan);
	}
	else
	{
		offsets = NULL;
	}

	if (it.scales != NULL)
	{
		if (scales != NULL) delete [] scales;
		scales = new float [npol * nchan];
		std::memcpy(scales, it.scales, sizeof(float) * npol * nchan);
	}
	else
	{
		scales = NULL;
	}

	if (it.data != NULL)
	{
		if (data != NULL) delete [] data;
		if (mode == FOLD)
		{
			data = new short [npol * nchan * nbin];
			std::memcpy(data, it.data, sizeof(short) * npol * nchan * nbin);
		}
		else if (mode == SEARCH)
		{
			switch (dtype)
			{
			case UINT1: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT2: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT4: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case UINT8: data = new unsigned char [nsblk * npol * nchan * nbits / 8]; std::memcpy(data, it.data, sizeof(unsigned char) * nsblk * npol * nchan * nbits / 8); break;
			case FLOAT: data = new float [nsblk * npol * nchan]; std::memcpy(data, it.data, sizeof(float) * nsblk * npol * nchan); break;
			default: cerr<<"Error: data type not support"<<endl; break;
			}
		}
	}
	else
	{
		data = NULL;
	}

	return *this;
}

Integration::~Integration()
{
	if (frequencies != NULL)
	{
		delete [] frequencies;
		frequencies = NULL;
	}

	if (weights != NULL)
	{
		delete [] weights;
		weights = NULL;
	}

	if (offsets != NULL)
	{
		delete [] offsets;
		offsets = NULL;
	}

	if (scales != NULL)
	{
		delete [] scales;
		scales = NULL;
	}

	if (data != NULL)
	{
		switch (dtype)
		{
		case UINT1: delete [] (unsigned char *)data; break;
		case UINT2: delete [] (unsigned char *)data; break;
		case UINT4: delete [] (unsigned char *)data; break;
		case UINT8: delete [] (unsigned char *)data; break;
		case FLOAT: delete [] (float *)data; break;
		case SHORT: delete [] (short *)data; break;
		case USHORT: delete [] (short *)data; break;
		default: cerr<<"Error: data type not support"<<endl; break;
		}
		data = NULL;
	}
}

void Integration::free()
{
	if (frequencies != NULL)
	{
		delete [] frequencies;
		frequencies = NULL;
	}

	if (weights != NULL)
	{
		delete [] weights;
		weights = NULL;
	}

	if (offsets != NULL)
	{
		delete [] offsets;
		offsets = NULL;
	}

	if (scales != NULL)
	{
		delete [] scales;
		scales = NULL;
	}

	if (data != NULL)
	{
		switch (dtype)
		{
		case UINT1: delete [] (unsigned char *)data; break;
		case UINT2: delete [] (unsigned char *)data; break;
		case UINT4: delete [] (unsigned char *)data; break;
		case UINT8: delete [] (unsigned char *)data; break;
		case FLOAT: delete [] (float *)data; break;
		case SHORT: delete [] (short *)data; break;
		case USHORT: delete [] (short *)data; break;
		default: cerr<<"Error: data type not support"<<endl; break;
		}
		data = NULL;
	}
}

void Integration::resize(int np, int nc, int nb)
{
	if (mode == FOLD)
	{
		if (np != npol or nc != nchan or nb != nbin)
		{
			npol = np;
			nchan = nc;
			nbin = nb;

			if (frequencies != NULL) delete [] frequencies;
			if (weights != NULL) delete [] weights;
			if (offsets != NULL) delete [] offsets;
			if (scales != NULL) delete [] scales;
			if (data != NULL)
			{
				delete [] (short *)data;
			}
			frequencies = new double [nchan];
			weights = new float [nchan];
			offsets = new float [npol*nchan];
			scales = new float [npol*nchan];
			data = new short [npol*nchan*nbin];
		}

		memset(data, 0, sizeof(short)*npol*nchan*nbin);
	}
	else if (mode == SEARCH)
	{
		if (np != npol or nc != nchan or nb != nsblk)
		{
			npol = np;
			nchan = nc;
			nsblk = nb;

			if (frequencies != NULL) delete [] frequencies;
			if (weights != NULL) delete [] weights;
			if (offsets != NULL) delete [] offsets;
			if (scales != NULL) delete [] scales;
			if (data != NULL)
			{
				switch (dtype)
				{
				case UINT1: delete [] (unsigned char *)data; break;
				case UINT2: delete [] (unsigned char *)data; break;
				case UINT4: delete [] (unsigned char *)data; break;
				case UINT8: delete [] (unsigned char *)data; break;
				case FLOAT: delete [] (float *)data; break;
				default: cerr<<"Error: data type not support"<<endl; break;
				}
			}
			frequencies = new double [nchan];
			weights = new float [nchan];
			offsets = new float [npol*nchan];
			scales = new float [npol*nchan];
			switch (dtype)
			{
			case UINT1: nbits = 1; data = new unsigned char [nsblk*npol*nchan*nbits/8];break;
			case UINT2: nbits = 2; data = new unsigned char [nsblk*npol*nchan*nbits/8]; break;
			case UINT4: nbits = 4; data = new unsigned char [nsblk*npol*nchan*nbits/8]; break;
			case UINT8: nbits = 8; data = new unsigned char [nsblk*npol*nchan*nbits/8]; break;
			case FLOAT: nbits = 32; data = new float [nsblk*npol*nchan]; break;
			default: cerr<<"Error: data type not support"<<endl; break;
			}
		}

		memset(data, 0, nsblk*npol*nchan*nbits/8);
	}

	for (long int i=0; i<npol*nchan; i++)
	{
		scales[i] = 1.;
		offsets[i] = 0.;
	}
	for (long int i=0; i<nchan; i++)
	{
		frequencies[i] = 0.;
		weights[i] = 1.;
	}
}

void Integration::load_data(void *dat, int np, int nc, int nb)
{
	resize(np, nc, nb);

	if (mode == FOLD)
	{
		long int TYPE_MAX = 0;
		long int TYPE_MIN = 0;
		if (dtype == SHORT)
		{
			TYPE_MAX = 32767;
			TYPE_MIN = 0;
		}
		else if (dtype == USHORT)
		{
			TYPE_MAX = 65535;
			TYPE_MIN = 0;
		}

		float *pro = (float *)dat;
		long int m = 0;
		long int n = 0;
		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				float max_value = pro[0];
				float min_value = pro[0];
				for (long int i=0; i<nbin; i++)
				{
					max_value = pro[i]>max_value ? pro[i]:max_value;
					min_value = pro[i]<min_value ? pro[i]:min_value;
				}
				scales[m] = (max_value-min_value)/(TYPE_MAX-TYPE_MIN);
				if (dtype == SHORT)
					offsets[m] = (max_value+min_value)*0.5;
				else if (dtype == USHORT)
					offsets[m] = min_value;
				for (long int i=0; i<nbin; i++)
				{
					((short *)data)[n] = (pro[i]-offsets[m])/scales[m];
					n++;
				}
				pro += nbin;
				m++;
			}
		}
		pro = NULL;
	}
	else if (mode == SEARCH)
	{
		if (dtype ==FLOAT)
		{
			memcpy(data, dat, sizeof(float)*nsblk*npol*nchan);
		}
		else
		{
			memcpy(data, dat, sizeof(unsigned char)*nsblk*npol*nchan*nbits/8);
		}
	}
}

void Integration::load_frequencies(double *freq, int nc)
{
	int nch = nc<nchan ? nc:nchan;

	for (long int j=0; j<nch; j++)
	{
		frequencies[j] = freq[j];
	}
}

void Integration::load_weights(float *wts, int nc)
{
	int nch = nc<nchan ? nc:nchan;

	for (long int j=0; j<nch; j++)
	{
		weights[j] = wts[j];
	}
}

bool Integration::to_char(Integration &it)
{
	assert(it.mode == SEARCH);
	assert(it.dtype == UINT8);
	assert(it.nbits == 8);
	assert(it.nsblk == nsblk);
	assert(it.npol == npol);
	assert(it.nchan == nchan);
	assert(it.data != NULL);

	switch (nbits)
	{
	case 8:
	{
		long int nchr = nsblk*npol*nchan;
		for (long int i=0; i<nchr; i++)
		{
			((unsigned char *)(it.data))[i] = ((unsigned char *)data)[i];
		}
	}
	break;

	case 4:
	{
		long int nchr = nsblk*npol*nchan;
		for (long int i=0; i<nchr/2; i++)
		{
			unsigned char tmp = ((unsigned char *)data)[i];
			for (long int k=0; k<2; k++)
			{
				((unsigned char *)(it.data))[i*2+k] = (tmp & 0b1111);
				tmp >>= 4;
			}
		}
	}
	break;
	
	case 2:
	{
		long int nchr = nsblk*npol*nchan;
		for (long int i=0; i<nchr/4; i++)
		{
			unsigned char tmp = ((unsigned char *)data)[i];
			for (long int k=0; k<4; k++)
			{
				((unsigned char *)(it.data))[i*4+k] = (tmp & 0b11);
				tmp >>= 2;
			}
		}
	}
	break;

	case 1:
	{
		long int nchr = nsblk*npol*nchan;
		for (long int i=0; i<nchr/8; i++)
		{
			unsigned char tmp = ((unsigned char *)data)[i];
			for (long int k=0; k<8; k++)
			{
				((unsigned char *)(it.data))[i*8+k] = tmp & 1;
				tmp >>= 1;
			}
		}
	}
	break;
	
	default:
	{
		std::cerr<<"Warning: data type unsupported"<<endl;
		return false;
	}
	break;
	}

	return true;
}