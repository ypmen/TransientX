/*
 * utils.h
 *
 *  Created on: Feb 26, 2020
 *      Author: ypmen
 */

#ifndef UTILS_H
#define UTILS_H

#include "config.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <set>
#ifdef __AVX2__
	#include <immintrin.h>
#endif
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include <random>
#include <boost/algorithm/string.hpp>

using namespace std;

template <typename T>
inline vector<size_t> argsort(const vector<T> &v)
{
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T>
inline vector<size_t> argsort2(const vector<T> &v)
{
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

template <typename T>
inline vector<size_t> argsort(const vector<vector<T>> &points, int dim)
{
  vector<size_t> idx(points.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(), [&points, dim](size_t p1, size_t p2) {return points[p1][dim] < points[p2][dim];});

  return idx;
}

long double to_longdouble(double value1, double value2);

template <typename T>
inline T kadane (T *arr, long int n, long int *start, long int *end)
{  
	T sum = 0;
	T maxSum = -std::numeric_limits<T>::infinity();
	*end = -1;

	long int local_start = 0;  

	for (long int i = 0; i < n; ++i)  
	{
		sum += arr[i];
		if (sum < 0)
		{
			sum = 0;
			local_start = i + 1;
		}
		else if (sum > maxSum)
		{
			maxSum = sum;
			*start = local_start;
			*end = i;
		}
	}

	if (*end != -1)
		return maxSum;

	maxSum = arr[0];
	*start = *end = 0;
  
	return maxSum;  
}

int gcd(int a, int b);
// Returns LCM of array elements 
long int findlcm(int arr[], int n);

fftwf_plan plan_transpose(int rows, int cols, float *in, float *out);
fftw_plan plan_transpose(int rows, int cols, double *in, double *out);

void runMedian(float *data, float *datMedian, long int size, int w);
template <typename T>
void runMedian2(T *data, T *datMedian, long int size, int w);
template <typename T>
void runMedian3(T *data, T *datMedian, long int size, int w);

#ifdef __AVX2__
void transpose_AVX2(float *out, float *in, int m, int n);
#endif

template <typename T>
void transpose(T *out, T *in, int m, int n);
template <typename T>
void transpose_pad(T *out, T *in, int m, int n);
template <typename T>
void transpose_pad(T *out, T *in, int m, int n, int tiley, int tilex);

void cmul(vector<complex<float>> &x, vector<complex<float>> &y);

template <typename T>
inline void dump2bin(const string &fname, const vector<T> &data)
{
	ofstream outfile;
	outfile.open(fname, ios::binary|ios::app);

	outfile.write((char *)(&data[0]), sizeof(T)*data.size());

	outfile.close();
}

/**
 * @brief Get the string ra and dec from double ra and dec
 * 
 * @param ra
 * @param dec
 * @param s_ra
 * @param s_dec
 */
void get_s_radec(double ra, double dec, string &s_ra, string &s_dec);

bool inverse_matrix4x4(const double m[16], double invOut[16]);
/**
 * @brief Get the inverse 4x4 matrix
 * 
 * @param m: input matrix
 * @param invOut: output inverse matrix
 * @return fasle if singluar 
 */
inline bool get_inverse_matrix4x4(const double m[4][4], double invOut[4][4])
{
	return inverse_matrix4x4((double *)m, (double *)invOut);
}

/**
 * @brief Get the inverse 3x3 matrix
 * 
 * @param m 
 * @param invOut
 * @return fasle if singluar 
 */
bool get_inverse_matrix3x3(const double m[3][3], double invOut[3][3]);

/**
 * @brief Get the error from chisq matrix
 * 
 * @param xerr: error of x
 * @param yerr: error of y
 * @param x:    shape = (N,)
 * @param y:    shape = (M,)
 * @param mxchisq:  chi square distribution with shape = (M, N)
 */
template <typename T>
bool get_error_from_chisq_matrix(T &xerr, T &yerr, vector<T> &x, vector<T> &y, vector<T> &mxchisq);
template <typename T>
bool get_error_from_chisq_matrix(T &xerr, vector<T> &x, vector<T> &vchisq);

/**
 * @brief format the output of value(err), such as 1.234(5), 1.23(4)e-6
 * 
 * @tparam T 
 * @param s_val_err 
 * @param val 
 * @param err 
 * @param style : sci or plain
 * @param low : number lower than low will be formatted to scientic, when style="sci"
 * @param high 
 */
template <typename T>
inline void format_val_err(std::string &s_val_err, T val, T err, const std::string &style="plain", int low=-5, int high=6)
{
	std::stringstream ss_val;
	std::stringstream ss_err;

	if (style == "sci" and (val<pow(10, low) or val>pow(10, high)))
	{
		int n = 0;
		if (val != 0)
		{
			n = floor(log10(abs(val)));
			val *= pow(10, -n);
			err *= pow(10, -n);
			
			if (err>=1)
			{
				ss_val<<fixed<<setprecision(0)<<val;
				ss_err<<fixed<<setprecision(0)<<err;
			}
			else
			{
				int n = -floor(log10(err))+1;
				ss_val<<fixed<<setprecision(n)<<val;
				ss_err<<fixed<<setprecision(0)<<err*pow(10, n);
			}
		}
		else
		{
			n = floor(log10(abs(err)));
			val *= pow(10, -n);
			err *= pow(10, -n);

			ss_val<<fixed<<setprecision(0)<<val;
			ss_err<<fixed<<setprecision(0)<<err;
		}
		
		s_val_err = ss_val.str() + "(" + ss_err.str() + ")" + "e"+to_string(n);
	}
	else
	{
		if (err>=1 or (val == 0 and err == 0))
		{
			ss_val<<fixed<<setprecision(0)<<val;
			ss_err<<fixed<<setprecision(0)<<err;
		}
		else
		{
			if ((err < abs(val) and err != 0) or val==0)
			{
				int n = -floor(log10(err))+1;
				ss_val<<fixed<<setprecision(n)<<val;
				ss_err<<fixed<<setprecision(0)<<err*pow(10, n);
			}
			else
			{
				int n = -floor(log10(abs(val)))+1;
				ss_val<<fixed<<setprecision(n)<<val;
				ss_err<<fixed<<setprecision(0)<<err*pow(10, n);
			}
		}

		s_val_err = ss_val.str() + "(" + ss_err.str() + ")";
	}    
}

/**
 * @brief Calculate ra dec in rad from string
 * 
 * @param s_ra 
 * @param s_dec 
 * @param ra 
 * @param dec 
 */
inline void get_rad_radec(const std::string &s_ra, const std::string &s_dec, double &ra, double &dec)
{
	std::vector<std::string> hhmmss;
	boost::split(hhmmss, s_ra, boost::is_any_of(":"), boost::token_compress_on);

	std::vector<std::string> ddmmss;
	boost::split(ddmmss, s_dec, boost::is_any_of(":"), boost::token_compress_on);

	for (long int i=hhmmss.size(); i<3; i++)
	{
		hhmmss.push_back("0");
	}
	for (long int i=ddmmss.size(); i<3; i++)
	{
		ddmmss.push_back("0");
	}

	double sign = std::signbit(stod(ddmmss[0])) ?  -1 : 1;
	ra = (stod(hhmmss[0]) + stod(hhmmss[1])/60. + stod(hhmmss[2])/3600.)*15./180.*M_PI;
	dec = sign*(sign*stod(ddmmss[0]) + stod(ddmmss[1])/60. + stod(ddmmss[2])/3600.)/180.*M_PI;
}

template <typename T>
T randnorm(const T &mean, const T &stddev)
{
	static thread_local std::mt19937 generator;
	std::normal_distribution<T> distribution(mean, stddev);
	return distribution(generator);
}

void get_gl_gb(double &gl, double &gb, const std::string &s_ra, const std::string &s_dec);
double get_maxdm_ymw16(double gl, double gb);
double get_dist_ymw16(double gl, double gb, double dm);

template <typename T>
void get_mean_var(T profile, int size, double &mean, double &var);

template <typename T>
void get_skewness_kurtosis(T profile, int size, double &skewness, double &kurtosis);

template <typename T>
void get_mean_var2(T profile, int size, double &mean, double &var);

template <typename T>
void get_mean_var(T profiles, int nrow, int ncol, double &mean, double &var);

template <typename T>
void get_mean_var(T profile, T profiles, int nsubint, int nchan, int nbin, double &mean, double &var);

inline void get_bestfit(float &a, float &b, const std::vector<float> &data, const std::vector<float> &data_ref)
{
	assert(data.size() == data_ref.size());
	int N = data.size();

	double xe = 0.;
	double ss = 0.;
	double ee = N;
	double se = 0.;
	double xs = 0.;

	for (long int i=0; i<N; i++)
	{
		xe += data[i];
		se += data_ref[i];
		ss += data_ref[i]*data_ref[i];
		xs += data[i]*data_ref[i];
	}

	double tmp = se*se-ss*ee;
	a = (xe*se-xs*ee)/tmp;
	b = (xs*se-xe*ss)/tmp;
}

#endif /* UTILS_H */
