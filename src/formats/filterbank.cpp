/*
 * filterbank.cpp
 *
 *  Created on: Feb 19, 2020
 *      Author: ypmen
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "filterbank.h"

using namespace std;

#define _CHAR_SWAP_SIZE 256

Filterbank::Filterbank()
{
	header_size = 0;
	use_frequence_table = false;

	telescope_id = 0;
	machine_id = 0;
	data_type = 1;
	barycentric = 0;
	pulsarcentric = 0;
	ibeam = 0;
	nbeams = 0;
	npuls = 0;
	nbins = 0;
	az_start = 0.;
	za_start = 0.;
	src_raj = 0.;
	src_dej = 0.;
	tstart = 0.;
	tsamp = 0.;
	nbits = 0;
	nsamples = 0.;
	nifs = 0;
	nchans = 0;
	fch1 = 0.;
	foff = 0.;
	refdm = 0.;
	period = 0.;

	frequency_table = new double [16320];
	ndata = 0;
	data = NULL;
	fptr = NULL;
}

Filterbank::Filterbank(const string fname)
{
	filename = fname;
	header_size = 0;
	use_frequence_table = false;

	telescope_id = 0;
	machine_id = 0;
	data_type = 1;
	barycentric = 0;
	pulsarcentric = 0;
	ibeam = 0;
	nbeams = 0;
	npuls = 0;
	nbins = 0;
	az_start = 0.;
	za_start = 0.;
	src_raj = 0.;
	src_dej = 0.;
	tstart = 0.;
	tsamp = 0.;
	nbits = 0;
	nsamples = 0.;
	nifs = 0;
	nchans = 0;
	fch1 = 0.;
	foff = 0.;
	refdm = 0.;
	period = 0.;

	frequency_table = new double [16320];
	ndata = 0;
	data = NULL;
	fptr = NULL;
}

Filterbank::Filterbank(const Filterbank &fil)
{
	filename = fil.filename;
	header_size = fil.header_size;
	use_frequence_table = fil.use_frequence_table;
	telescope_id = fil.telescope_id;
	machine_id = fil.machine_id;
	data_type = fil.data_type;
	strcpy(rawdatafile, fil.rawdatafile);
	strcpy(source_name, fil.source_name);
	barycentric = fil.barycentric;
	pulsarcentric = fil.pulsarcentric;
	ibeam = fil.ibeam;
	nbeams = fil.nbeams;
	npuls = fil.npuls;
	nbins = fil.nbins;
	az_start = fil.az_start;
	za_start = fil.za_start;
	src_raj = fil.src_raj;
	src_dej = fil.src_dej;
	tstart = fil.tstart;
	tsamp = fil.tsamp;
	nbits = fil.nbits;
	nsamples = fil.nsamples;
	nifs = fil.nifs;
	nchans = fil.nchans;
	fch1 = fil.fch1;
	foff = fil.foff;
	refdm = fil.refdm;
	period = fil.period;

	if (fil.frequency_table != NULL)
	{
		frequency_table = new double [16320];
		memcpy(frequency_table, fil.frequency_table, sizeof(double)*16320);
	}
	else
	{
		frequency_table = NULL;
	}

	ndata = fil.ndata;

	if (fil.data != NULL)
	{
		switch (nbits)
		{
		case 1:
		{
			data = new unsigned char [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(unsigned char)*ndata*nifs*nchans); break;
		}
		case 8:
		{
			data = new unsigned char [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(unsigned char)*ndata*nifs*nchans); break;
		}
		case 32:
		{
			data = new float [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(float)*ndata*nifs*nchans); break;
		}
		default: cerr<<"Error: data type not support"<<endl; break;
		}
	}
	else
	{
		data = NULL;
	}

	fptr = NULL;
}

Filterbank & Filterbank::operator=(const Filterbank &fil)
{
	filename = fil.filename;
	header_size = fil.header_size;
	use_frequence_table = fil.use_frequence_table;
	telescope_id = fil.telescope_id;
	machine_id = fil.machine_id;
	data_type = fil.data_type;
	strcpy(rawdatafile, fil.rawdatafile);
	strcpy(source_name, fil.source_name);
	barycentric = fil.barycentric;
	pulsarcentric = fil.pulsarcentric;
	ibeam = fil.ibeam;
	nbeams = fil.nbeams;
	npuls = fil.npuls;
	nbins = fil.nbins;
	az_start = fil.az_start;
	za_start = fil.za_start;
	src_raj = fil.src_raj;
	src_dej = fil.src_dej;
	tstart = fil.tstart;
	tsamp = fil.tsamp;
	nbits = fil.nbits;
	nsamples = fil.nsamples;
	nifs = fil.nifs;
	nchans = fil.nchans;
	fch1 = fil.fch1;
	foff = fil.foff;
	refdm = fil.refdm;
	period = fil.period;

	if (fil.frequency_table != NULL)
	{
		if (frequency_table != NULL) delete [] frequency_table;
		frequency_table = new double [16320];
		memcpy(frequency_table, fil.frequency_table, sizeof(double)*16320);
	}

	ndata = fil.ndata;

	if (fil.data != NULL)
	{
		switch (nbits)
		{
		case 1:
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(unsigned char)*ndata*nifs*nchans); break;
		}
		case 8:
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(unsigned char)*ndata*nifs*nchans); break;
		}
		case 32:
		{
			if (data != NULL) delete [] (float *)data;
			data = new float [ndata*nifs*nchans];
			memcpy(data, fil.data, sizeof(float)*ndata*nifs*nchans); break;
		}
		default: cerr<<"Error: data type not support"<<endl; break;
		}
	}

	return *this;
}

Filterbank::~Filterbank()
{
	if (frequency_table != NULL)
	{
		delete [] frequency_table;
		frequency_table = NULL;
	}

	if (data != NULL)
	{
		switch (nbits)
		{
		case 1: delete [] (unsigned char *)data; break;
		case 2: delete [] (unsigned char *)data; break;
		case 4: delete [] (unsigned char *)data; break;
		case 8: delete [] (unsigned char *)data; break;
		case 32: delete [] (float *)data; break;
		default: cerr<<"Error: data type not support"<<endl; break;
		}
		data = NULL;
	}

	if (fptr != NULL)
	{
		fclose(fptr);
		fptr = NULL;
	}
}

void Filterbank::free()
{
	if (frequency_table != NULL)
	{
		delete [] frequency_table;
		frequency_table = NULL;
	}

	if (data != NULL)
	{
		switch (nbits)
		{
		case 1: delete [] (unsigned char *)data; break;
		case 2: delete [] (unsigned char *)data; break;
		case 4: delete [] (unsigned char *)data; break;
		case 8: delete [] (unsigned char *)data; break;
		case 32: delete [] (float *)data; break;
		default: cerr<<"Error: data type not support"<<endl; break;
		}
		data = NULL;
	}

	if (fptr != NULL)
	{
		fclose(fptr);
		fptr = NULL;
	}
}

void Filterbank::close()
{
	if (fptr != NULL)
	{
		fclose(fptr);
		fptr = NULL;
	}
}

bool Filterbank::read_header()
{
	int expecting_rawdatafile = 0, expecting_source_name = 0;
	int nsamp;
	int expecting_frequency_table = 0, channel_index = 0;
	fptr = fopen(filename.c_str (), "rb");
	if (!fptr)
	{
		cerr<<"Can not open file."<<endl;
		return false;
	}

	string strtmp;
	long int intTotalHeaderBytes = 0;
	get_string(fptr, strtmp);
	if (strtmp != "HEADER_START")
	{
		cerr << "Non-Standard file format." << endl;
		return false;
	}
	intTotalHeaderBytes += strtmp.length () + 4;
	while (1)
	{
		get_string (fptr, strtmp);
		if (strtmp == "HEADER_END")
		{
			intTotalHeaderBytes += strtmp.length () + 4;
			break;
		}
		intTotalHeaderBytes += strtmp.length () + 4;
		if (strtmp == "rawdatafile")
		{
			expecting_rawdatafile = 1;
		}
		else if (strtmp == "source_name")
		{
			expecting_source_name = 1;
		}
		else if (strtmp == "FREQUENCY_START")
		{
			expecting_frequency_table = 1;
			use_frequence_table=true;
			channel_index = 0;
		}
		else if (strtmp == "FREQUENCY_END")
		{
			expecting_frequency_table = 0;
			nchans=channel_index;
		}
		else if (strtmp == "az_start")
		{
			fread(&az_start, sizeof(az_start), 1, fptr);
			intTotalHeaderBytes += sizeof(az_start);
		}
		else if (strtmp == "za_start")
		{
			fread(&za_start, sizeof(za_start), 1, fptr);
			intTotalHeaderBytes += sizeof(za_start);
		}
		else if (strtmp == "src_raj")
		{
			fread(&src_raj, sizeof(src_raj), 1, fptr);
			intTotalHeaderBytes += sizeof(src_raj);
		}
		else if (strtmp == "src_dej")
		{
			fread(&src_dej, sizeof(src_dej), 1, fptr);
			intTotalHeaderBytes += sizeof(src_dej);
		}
		else if (strtmp == "tstart")
		{
			fread(&tstart, sizeof(tstart), 1, fptr);
			intTotalHeaderBytes += sizeof(tstart);
		}
		else if (strtmp == "tsamp")
		{
			fread(&tsamp, sizeof(tsamp), 1, fptr);
			intTotalHeaderBytes += sizeof(tsamp);
		}
		else if (strtmp == "period")
		{
			fread(&period, sizeof(period), 1, fptr);
			intTotalHeaderBytes += sizeof(period);
		}
		else if (strtmp == "fch1")
		{
			fread(&fch1, sizeof(fch1), 1, fptr);
			intTotalHeaderBytes += sizeof(fch1);
		}
		else if (strtmp == "fchannel")
		{
			fread(&frequency_table[channel_index], sizeof(double), 1, fptr);
			intTotalHeaderBytes += sizeof(double);
			fch1 = foff = 0.0;				/* set to 0.0 to signify that a table is in use */
			use_frequence_table = true;
			channel_index++;
		}
		else if (strtmp == "foff")
		{
			fread(&foff, sizeof(foff), 1, fptr);
			intTotalHeaderBytes += sizeof(foff);
		}
		else if (strtmp == "nchans")
		{
			fread(&nchans, sizeof(nchans), 1, fptr);
			intTotalHeaderBytes += sizeof(nchans);
		}
		else if (strtmp == "telescope_id")
		{
			fread(&telescope_id, sizeof(telescope_id), 1, fptr);
			intTotalHeaderBytes += sizeof(telescope_id);
		}
		else if (strtmp == "machine_id")
		{
			fread(&machine_id, sizeof(machine_id), 1, fptr);
			intTotalHeaderBytes += sizeof(machine_id);
		}
		else if (strtmp == "data_type")
		{
			fread(&data_type, sizeof(data_type), 1, fptr);
			intTotalHeaderBytes += sizeof(data_type);
		}
		else if (strtmp == "ibeam")
		{
			fread(&ibeam, sizeof(ibeam), 1, fptr);
			intTotalHeaderBytes += sizeof(ibeam);
		}
		else if (strtmp == "nbeams")
		{
			fread(&nbeams, sizeof(nbeams), 1, fptr);
			intTotalHeaderBytes += sizeof(nbeams);
		}
		else if (strtmp == "nbits")
		{
			fread(&nbits, sizeof(nbits), 1, fptr);
			intTotalHeaderBytes += sizeof(nbits);
		}
		else if (strtmp == "barycentric")
		{
			fread(&barycentric, sizeof(barycentric), 1, fptr);
			intTotalHeaderBytes += sizeof(barycentric);
		}
		else if (strtmp == "pulsarcentric")
		{
			fread(&pulsarcentric, sizeof(pulsarcentric), 1, fptr);
			intTotalHeaderBytes += sizeof(pulsarcentric);
		}
		else if (strtmp == "nbins")
		{
			fread(&nbins, sizeof(nbins), 1, fptr);
			intTotalHeaderBytes += sizeof(nbins);
		}
		else if (strtmp == "nsamples")
		{
			/* read this one only for backwards compatibility */
			fread(&nsamp, sizeof(nsamp), 1, fptr);
			intTotalHeaderBytes += sizeof(nsamp);
		}
		else if (strtmp == "nifs")
		{
			fread(&nifs, sizeof(nifs), 1, fptr);
			intTotalHeaderBytes += sizeof(nifs);
		}
		else if (strtmp == "npuls")
		{
			fread(&npuls, sizeof(npuls), 1, fptr);
			intTotalHeaderBytes += sizeof(npuls);
		}
		else if (strtmp == "refdm")
		{
			fread(&refdm, sizeof(refdm), 1, fptr);
			intTotalHeaderBytes += sizeof(refdm);
		}
		else if (expecting_rawdatafile == 1)
		{
			strcpy(rawdatafile, strtmp.c_str());
			expecting_rawdatafile = 0;
		}
		else if (expecting_source_name == 1)
		{
			strcpy(source_name, strtmp.c_str());
			expecting_source_name = 0;
		}
		else
		{
			cerr << "read_header - unknown parameter : " << strtmp << endl;
			return (false);
		}
	}
	header_size=intTotalHeaderBytes;
	nsamples=get_nsamples(filename.c_str(), header_size, nbits, nifs, nchans);
	//Correct Frequency table any way
	if (!use_frequence_table)
	{
		for (long int i=0; i<nchans; i++)
		{
			frequency_table[i]=fch1+i*foff;
		}
	}

	return true;
}

bool Filterbank::read_data()
{
	switch (nbits)
	{
	case 8:
	{
		long int nchr = nsamples*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr, fptr);
		data = new unsigned char [icnt];
		for (long int i=0; i<icnt; i++)
		{
				((unsigned char *)data)[i] = chb[i];
		}
		delete [] chb;
		if (icnt != nchr)
		{
				cerr<<"Data ends unexpected read to EOF"<<endl;
		}
		nsamples = icnt/nifs/nchans;
		ndata = nsamples;
	}; break;
	case 1:
	{
		long int nchr = nsamples*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/8, fptr);
		if (icnt*8 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*8];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<8; k++)
			{
				((unsigned char *)data)[i*8+k] = tmp & 1;
				tmp >>= 1;
			}
		}
		delete [] chb;
		if (icnt*8 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*8/nifs/nchans;
	}; break;
	case 2:
	{
		long int nchr = nsamples*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/4, fptr);
		if (icnt*4 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*4];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<4; k++)
			{
				((unsigned char *)data)[i*4+k] = (tmp & 0b11);
				tmp >>= 2;
			}
		}
		delete [] chb;
		if (icnt*4 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*4/nifs/nchans;
	}; break;
	case 4:
	{
		long int nchr = nsamples*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/2, fptr);
		if (icnt*2 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*2];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<2; k++)
			{
				((unsigned char *)data)[i*2+k] = (tmp & 0b1111);
				tmp >>= 4;
			}
		}
		delete [] chb;
		if (icnt*2 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*2/nifs/nchans;
	}; break;
	default:
	{
		cerr<<"Error: data type unsupported"<<endl;
		return false;
	}; break;
	}

	return true;
}

bool Filterbank::read_data(long int nstart, long int ns)
{
	long int offset=(long int) (long double)nstart * (long double)nchans * (long double)nifs * ((long double)nbits / 8.0);
	fseek(fptr, offset+header_size, SEEK_SET);
	return read_data(ns);
}

bool Filterbank::read_data(long int ns)
{
	switch (nbits)
	{
	case 8:
	{
		long int nchr = ns*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr, fptr);
		if (icnt > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt];
		}
		for (long int i=0; i<icnt; i++)
		{
				((unsigned char *)data)[i] = chb[i];
		}
		delete [] chb;
		if (icnt != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt/nifs/nchans;
	}; break;
	case 1:
	{
		long int nchr = ns*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/8, fptr);
		if (icnt*8 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*8];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<8; k++)
			{
				((unsigned char *)data)[i*8+k] = tmp & 1;
				tmp >>= 1;
			}
		}
		delete [] chb;
		if (icnt*8 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*8/nifs/nchans;
	}; break;
	case 2:
	{
		long int nchr = ns*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/4, fptr);
		if (icnt*4 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*4];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<4; k++)
			{
				((unsigned char *)data)[i*4+k] = (tmp & 0b11);
				tmp >>= 2;
			}
		}
		delete [] chb;
		if (icnt*4 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*4/nifs/nchans;
	}; break;
	case 4:
	{
		long int nchr = ns*nifs*nchans;
		unsigned char * chb = new unsigned char [nchr];
		long int icnt = fread(chb, 1, nchr/2, fptr);
		if (icnt*2 > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [icnt*2];
		}
		for (long int i=0; i<icnt; i++)
		{
			unsigned char tmp = chb[i];
			for (long int k=0; k<2; k++)
			{
				((unsigned char *)data)[i*2+k] = (tmp & 0b1111);
				tmp >>= 4;
			}
		}
		delete [] chb;
		if (icnt*2 != nchr)
		{
				//cerr<<"Warning: Data ends unexpected read to EOF"<<endl;
		}
		ndata = icnt*2/nifs/nchans;
	}; break;
	default:
	{
		cerr<<"Error: data type unsupported"<<endl;
		return false;
	}; break;
	}

	return true;
}

bool Filterbank::set_data(unsigned char *dat, long int ns, int nif, int nchan)
{
	switch (nbits)
	{
	case 8:
	{
		nifs = nif;
		nchans = nchan;
		long int nchr = ns*nifs*nchans*nbits/8;
		if (ns > ndata)
		{
			if (data != NULL) delete [] (unsigned char *)data;
			data = new unsigned char [nchr];
		}
		for (long int i=0; i<nchr; i++)
		{
				((unsigned char *)data)[i] = dat[i];
		}
		ndata = ns;
	}; break;
	default:
	{
		cerr<<"Error: data type unsupported"<<endl;
		return false;
	}; break;
	}

	return true;
}

bool Filterbank::write_header()
{
	fptr = fopen(filename.c_str(), "wb");
	if (fptr==NULL) return false;

	put_string(fptr, "HEADER_START");
	put_string(fptr, "source_name");
	put_string(fptr, source_name);
	if (use_frequence_table || (fch1==0.0 && foff==0.0))
	{
			put_string(fptr,"FREQUENCY_START");
			for (int channel_index=0; channel_index<nchans; channel_index++)
			{
					put_string(fptr, "fchannel");
					fwrite (&frequency_table[channel_index], sizeof(double), 1, fptr);
			}
			put_string(fptr, "FREQUENCY_END");
	}
	put_string(fptr, "az_start");
	fwrite(&az_start, sizeof(az_start), 1, fptr);
	put_string(fptr, "za_start");
	fwrite(&za_start, sizeof(za_start), 1, fptr);
	put_string(fptr, "src_raj");
	fwrite(&src_raj, sizeof(src_raj), 1, fptr);
	put_string(fptr, "src_dej");
	fwrite(&src_dej, sizeof(src_dej), 1, fptr);
	put_string(fptr, "tstart");
	fwrite(&tstart, sizeof(tstart), 1, fptr);
	put_string(fptr, "tsamp");
	fwrite(&tsamp, sizeof(tsamp), 1, fptr);
	//put_string(fptr, "period");
	//fwrite (&period, sizeof(period), 1, fptr);
	put_string(fptr, "fch1");
	fwrite(&fch1, sizeof(fch1), 1, fptr);
	put_string(fptr, "foff");
	fwrite(&foff, sizeof(foff), 1, fptr);
	put_string(fptr, "nchans");
	fwrite(&nchans, sizeof(nchans), 1, fptr);
	put_string(fptr, "telescope_id");
	fwrite(&telescope_id, sizeof(telescope_id), 1, fptr);
	put_string(fptr, "machine_id");
	fwrite(&machine_id, sizeof(machine_id), 1, fptr);
	put_string(fptr, "data_type");
	fwrite(&data_type, sizeof(data_type), 1, fptr);
	put_string(fptr, "ibeam");
	fwrite(&ibeam, sizeof(ibeam), 1, fptr);
	put_string(fptr, "nbeams");
	fwrite(&nbeams, sizeof(nbeams), 1, fptr);
	put_string(fptr, "nbits");
	fwrite(&nbits, sizeof(nbits), 1, fptr);
	put_string(fptr, "barycentric");
	fwrite(&barycentric, sizeof(barycentric), 1, fptr);
	put_string(fptr, "pulsarcentric");
	fwrite(&pulsarcentric, sizeof(pulsarcentric), 1, fptr);
	//put_string(fptr, "nbins");
	//fwrite (&nbins, sizeof(nbins), 1, fptr);
	//put_string(fptr, "nsamples");
	//fwrite (&nsamples, sizeof(nsamples), 1, fptr);
	put_string(fptr, "nifs");
	fwrite(&nifs, sizeof(nifs), 1, fptr);
	//put_string(fptr, "npuls");
	//fwrite (&npuls, sizeof(npuls), 1, fptr);
	if (data_type == 2)
	{
		put_string(fptr, "refdm");
		fwrite (&refdm, sizeof(refdm), 1, fptr);
	}
	put_string(fptr, "HEADER_END");
	return true;
}

bool Filterbank::write_data()
{
	switch (nbits)
	{
	case 8:
	{
		long int nchr=ndata*nifs*nchans*nbits/8;
		fwrite(data, 1, nchr, fptr);
	}; break;
	case 32:
	{
		long int nchr=ndata*nifs*nchans;
		fwrite(data, 4, nchr, fptr);
	}; break;
	default: cerr<<"Error: data type is not supported"<<endl;
	}

	return true;
}

void Filterbank::put_string (FILE * outputfile, const string & strtmp)
{
		int nchar = strtmp.length ();
		fwrite (&nchar, sizeof (int), 1, outputfile);
		fwrite (strtmp.c_str (), nchar, 1, outputfile);
}

void Filterbank::get_string (FILE * inputfile, string & strtmp)
{
	int nchar;
	strtmp = "";
	fread (&nchar, sizeof (int), 1, inputfile);
	if (feof (inputfile))
	{
		cerr << "Error in reading the header of File" << endl;
		exit (0);
	}
	if (nchar > 80 || nchar < 1)
		return;
	char chrtmp[_CHAR_SWAP_SIZE];
	fread (chrtmp, nchar, 1, inputfile);
	chrtmp[nchar] = '\0';
	strtmp = chrtmp;
}

int Filterbank::get_nsamples(const char *filename,int headersize, int nbits, int nifs, int nchans)
{
	long long datasize;
	int numsamps;
	datasize=sizeof_file(filename)-headersize;
	numsamps=(int) ((long double) (datasize)/ (((long double) nbits) / 8.0)
					/(long double) nifs/(long double) nchans);
	return (numsamps);
}

long long Filterbank::sizeof_file(const char name[]) /* includefile */
{
	struct stat stbuf;

	if(stat(name,&stbuf) == -1)
	{
		fprintf(stderr, "f_siz: can't access %s\n",name);
		exit(0);
	}
	return(stbuf.st_size);
}

//copy from presto
void get_telescope_name(int telescope_id, std::string &s_telescope)
{
	switch (telescope_id) {
	case 0:
		s_telescope = "Fake";
		break;
	case 1:
		s_telescope = "Arecibo";
		break;
	case 2:
		s_telescope = "Ooty";
		break;
	case 3:
		s_telescope = "Nancay";
		break;
	case 4:
		s_telescope = "Parkes";
		break;
	case 5:
		s_telescope = "Jodrell";
		break;
	case 6:
		s_telescope = "GBT";
		break;
	case 7:
		s_telescope = "GMRT";
		break;
	case 8:
		s_telescope = "Effelsberg";
		break;
	case 9:
		s_telescope = "ATA";
		break;
	case 10:
		s_telescope = "SRT";
		break;
	case 11:
		s_telescope = "LOFAR";
		break;
	case 12:
		s_telescope = "VLA";
		break;
	case 20:  // May need to change....
		s_telescope = "CHIME";
		break;
	case 21:  // May need to change....
		s_telescope = "FAST";
		break;
	case 64:
		s_telescope = "MeerKAT";
		break;
	case 65:
		s_telescope = "KAT-7";
		break;
	default:
		s_telescope = "Unknown";
		break;
	}
}

bool iequals(const string& a, const string& b)
{
	unsigned int sz = a.size();
	if (b.size() != sz)
		return false;
	for (unsigned int i = 0; i < sz; ++i)
		if (tolower(a[i]) != tolower(b[i]))
			return false;
	return true;
}

int get_telescope_id(const std::string &s_telescope)
{
	if (iequals(s_telescope, "Fake"))
	{
		return 0;
	}
	else if (iequals(s_telescope, "Arecibo"))
	{
		return 1;
	}
	else if (iequals(s_telescope, "Ooty"))
	{
		return 2;
	}
	else if (iequals(s_telescope, "Nancay"))
	{
		return 3;
	}
	else if (iequals(s_telescope, "Parkes"))
	{
		return 4;
	}
	else if (iequals(s_telescope, "Jodrell"))
	{
		return 5;
	}
	else if (iequals(s_telescope, "GBT"))
	{
		return 6;
	}
	else if (iequals(s_telescope, "GMRT"))
	{
		return 7;
	}
	else if (iequals(s_telescope, "Effelsberg"))
	{
		return 8;
	}
	else if (iequals(s_telescope, "ATA"))
	{
		return 9;
	}
	else if (iequals(s_telescope, "SRT"))
	{
		return 10;
	}
	else if (iequals(s_telescope, "LOFAR"))
	{
		return 11;
	}
	else if (iequals(s_telescope, "VLA"))
	{
		return 12;
	}
	else if (iequals(s_telescope, "CHIME"))
	{
		return 20;
	}
	else if (iequals(s_telescope, "FAST"))
	{
		return 21;
	}
	else if (iequals(s_telescope, "MeerKAT"))
	{
		return 64;
	}
	else if (iequals(s_telescope, "KAT-7"))
	{
		return 65;
	}
	
	return -1;
}
