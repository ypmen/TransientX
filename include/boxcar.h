/*
 * boxcar.h
 *
 *  Created on: Apr 29, 2020
 *      Author: ypmen
 */

#ifndef BOXCAR_H_
#define BOXCAR_H_

#include <vector>
#include <utility>

#include "dedisperse.h"
#include "subdedispersion.h"

class Boxcar
{
public:
	Boxcar();
	Boxcar(const Boxcar &boxcar);
	Boxcar & operator=(const Boxcar &boxcar);
	~Boxcar();
	void prepare(RealTime::SubbandDedispersion &dedisp);
	void resize(long int nt, long int ndm);
	bool run(RealTime::SubbandDedispersion &dedisp, vector<int> &vwn, bool iqr=false);
	void match(int idm, vector<int> &vwn, RealTime::SubbandDedispersion &dedisp, bool iqr);
	void match2D(int idm, vector<int> &vwn, RealTime::SubbandDedispersion &dedisp);
	void dump2txt(const string fname) const;
public:
	long int counter;
	double tsamp;
	long int nsamples;
	double dms;
	double ddm;
	long int ndm;
	int maxwn;
	double fmin;
	double fmax;
	float *mxS;
	int *mxwn;
	vector<vector<std::pair<int, int>>> finterval;
};

#endif /* BOXCAR_H_ */
