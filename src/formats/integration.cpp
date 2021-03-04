/*
 * integration.cpp
 *
 *  Created on: Feb 22, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <cstring>

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
