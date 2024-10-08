/*
 * dedisperse.h
 *
 *  Created on: Apr 29, 2020
 *      Author: ypmen
 */

#ifndef DEDISPERSE_H_
#define DEDISPERSE_H_

#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils.h"
#include "json.hpp"

#define LOWLEVEL 0.1  //[0-1]
#define HIGHLEVEL 0.9 //[0-1]
#define HAVE_YMW16 1

extern unsigned int num_threads;

#endif /* DEDISPERSE_H_ */
