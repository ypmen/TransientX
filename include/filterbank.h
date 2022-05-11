/*
 * filterbank.h
 *
 *  Created on: Feb 19, 2020
 *      Author: ypmen
 */

#ifndef FILTERBANK_H_
#define FILTERBANK_H_

#include <stdio.h>
#include <string.h>

using namespace std;

class Filterbank
{
public:
	Filterbank();
	Filterbank(const Filterbank &fil);
	Filterbank & operator=(const Filterbank &fil);
	Filterbank(const string fname);
	~Filterbank();
	void free();
	void close();
	bool read_header();
	bool read_data();
	bool read_data(long int nstart, long int ns);
	bool read_data(long int ns);
	bool set_data(unsigned char *dat, long int ns, int nif, int nchan);
	bool write_header();
	bool write_data();
private:
	static void put_string (FILE * outputfile, const string & strtmp);
	static void get_string(FILE * inputfile, string & strtmp);
	static int get_nsamples(const char *filename, int headersize, int nbits, int nifs, int nchans);
	static long long sizeof_file(const char name[]);

public:
	string filename;
	long int header_size;
	bool use_frequence_table;

	int telescope_id;
	int machine_id;
	int data_type;
	char rawdatafile[80];
	char source_name[80];
	int barycentric;
	int pulsarcentric;
	int ibeam;
	int nbeams;
	int npuls;
	int nbins;
	double az_start;
	double za_start;
	double src_raj;
	double src_dej;
	double tstart;
	double tsamp;
	int nbits;
	long int nsamples;
	int nifs;
	int nchans;
	double fch1;
	double foff;
	double refdm;
	double period;

	double *frequency_table;
	long int ndata;
	void *data;
	FILE *fptr;
};

void get_telescope_name(int telescope_id, std::string &s_telescope);
int get_telescope_id(const std::string &s_telescope);

#endif /* FILTERBANK_H_ */
