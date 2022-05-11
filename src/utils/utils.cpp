/*
 * utils.cpp
 *
 *  Created on: Feb 27, 2020
 *      Author: ypmen
 */

#include <limits>
#include <assert.h>
#include <set>
#include "utils.h"
#include "AVL.h"
#include "dedisperse.h"

long double to_longdouble(double value1, double value2)
{
	return (long double)(value1) + (long double)(value2);
}

int gcd(int a, int b)
{
	if (b == 0)
		return a;
	return gcd(b, a % b);
}

// Returns LCM of array elements
long int findlcm(int arr[], int n)
{
	// Initialize result
	long int ans = arr[0];

	// ans contains LCM of arr[0], ..arr[i]
	// after i'th iteration,
	for (int i = 1; i < n; i++)
		ans = (((arr[i] * ans)) /
			   (gcd(arr[i], ans)));
	return ans;
}

fftwf_plan plan_transpose(int rows, int cols, float *in, float *out)
{
	const unsigned flags = FFTW_ESTIMATE; /* other flags are possible */
	fftw_iodim howmany_dims[2];

	howmany_dims[0].n = rows;
	howmany_dims[0].is = cols;
	howmany_dims[0].os = 1;

	howmany_dims[1].n = cols;
	howmany_dims[1].is = 1;
	howmany_dims[1].os = rows;

	return fftwf_plan_guru_r2r(/*rank=*/0, /*dims=*/NULL,
							   /*howmany_rank=*/2, howmany_dims,
							   in, out, /*kind=*/NULL, flags);
}

fftw_plan plan_transpose(int rows, int cols, double *in, double *out)
{
	const unsigned flags = FFTW_ESTIMATE; /* other flags are possible */
	fftw_iodim howmany_dims[2];

	howmany_dims[0].n = rows;
	howmany_dims[0].is = cols;
	howmany_dims[0].os = 1;

	howmany_dims[1].n = cols;
	howmany_dims[1].is = 1;
	howmany_dims[1].os = rows;

	return fftw_plan_guru_r2r(/*rank=*/0, /*dims=*/NULL,
							  /*howmany_rank=*/2, howmany_dims,
							  in, out, /*kind=*/NULL, flags);
}

template <typename T>
void transpose(T *out, T *in, int m, int n)
{
	const int tilex = 16;
	const int tiley = 64;

	assert(n % tilex == 0);
	assert(m % tiley == 0);

	int blockx = n / tilex;
	int blocky = m / tiley;

	T *temp = new T[num_threads * tiley * tilex];
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int s = 0; s < blocky*blockx; s++)
	{
		T *ptemp = temp;
#ifdef _OPENMP
		ptemp = temp + omp_get_thread_num()*tiley*tilex;
#endif
		long int l = s/blockx;
		long int k = s%blockx;

		for (long int i = 0; i < tiley; i++)
		{
			for (long int j = 0; j < tilex; j++)
			{
				ptemp[j * tiley + i] = in[l * tiley * n + i * n + k * tilex + j];
			}
		}
		for (long int i = 0; i < tilex; i++)
		{
			for (long int j = 0; j < tiley; j++)
			{
				out[k * tilex * m + i * m + l * tiley + j] = ptemp[i * tiley + j];
			}
		}
	}

	delete[] temp;
}

template <typename T>
void transpose_pad(T *out, T *in, int m, int n)
{
	const int tilex = 16;
	const int tiley = 64;

	int npad = ceil(n*1./tilex)*tilex;
	int mpad = ceil(m*1./tiley)*tiley;

	int blockx = npad / tilex;
	int blocky = mpad / tiley;

	T *temp = new T[num_threads * tiley * tilex];
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int s = 0; s < blocky*blockx; s++)
	{
		T *ptemp = temp;
#ifdef _OPENMP
		ptemp = temp + omp_get_thread_num()*tiley*tilex;
#endif
		long int l = s/blockx;
		long int k = s%blockx;
		
		if ((l+1)*tiley <= m and (k+1)*tilex <= n)
		{
			for (long int i = 0; i < tiley; i++)
			{
				for (long int j = 0; j < tilex; j++)
				{
					ptemp[j * tiley + i] = in[(l * tiley + i) * n + k * tilex + j];
				}
			}
			for (long int i = 0; i < tilex; i++)
			{
				for (long int j = 0; j < tiley; j++)
				{
					out[(k * tilex + i) * m + l * tiley + j] = ptemp[i * tiley + j];
				}
			}
		}
		else if ((l+1)*tiley > m and (k+1)*tilex <= n)
		{
			for (long int i = 0; i < m-l*tiley; i++)
			{
				for (long int j = 0; j < tilex; j++)
				{
					ptemp[j * tiley + i] = in[(l * tiley + i) * n + k * tilex + j];
				}
			}
			for (long int i = 0; i < tilex; i++)
			{
				for (long int j = 0; j < m-l*tiley; j++)
				{
					out[(k * tilex + i) * m + l * tiley + j] = ptemp[i * tiley + j];
				}
			}
		}
		else if ((l+1)*tiley <= m and (k+1)*tilex > n)
		{
			for (long int i = 0; i < tiley; i++)
			{
				for (long int j = 0; j < n-k*tilex; j++)
				{
					ptemp[j * tiley + i] = in[(l * tiley + i) * n + k * tilex + j];
				}
			}
			for (long int i = 0; i < n-k*tilex; i++)
			{
				for (long int j = 0; j < tiley; j++)
				{
					out[(k * tilex + i) * m + l * tiley + j] = ptemp[i * tiley + j];
				}
			}
		}
		else
		{
			for (long int i = 0; i < m-l*tiley; i++)
			{
				for (long int j = 0; j < n-k*tilex; j++)
				{
					ptemp[j * tiley + i] = in[(l * tiley + i) * n + k * tilex + j];
				}
			}
			for (long int i = 0; i < n-k*tilex; i++)
			{
				for (long int j = 0; j < m-l*tiley; j++)
				{
					out[(k * tilex + i) * m + l * tiley + j] = ptemp[i * tiley + j];
				}
			}
		}
	}

	delete[] temp;
}

template <typename T>
void transpose_pad(T *out, T *in, int m, int n, int tiley, int tilex)
{
	int npad = ceil(n*1./tilex)*tilex;
	int mpad = ceil(m*1./tiley)*tiley;

	int blockx = npad / tilex;
	int blocky = mpad / tiley;

	T *temp = new T[num_threads * tiley * tilex];
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int s = 0; s < blocky*blockx; s++)
	{
		T *ptemp = temp;
#ifdef _OPENMP
		ptemp = temp + omp_get_thread_num()*tiley*tilex;
#endif
		long int l = s/blockx;
		long int k = s%blockx;
		
		for (long int i = 0; i < tiley; i++)
		{
			for (long int j = 0; j < tilex; j++)
			{
				if (l * tiley + i < m and k * tilex + j < n)
					ptemp[j * tiley + i] = in[(l * tiley + i) * n + k * tilex + j];
			}
		}
		for (long int i = 0; i < tilex; i++)
		{
			for (long int j = 0; j < tiley; j++)
			{
				if (k * tilex + i < n and l * tiley + j < m)
					out[(k * tilex + i) * m + l * tiley + j] = ptemp[i * tiley + j];
			}
		}
	}

	delete[] temp;    
}

#ifdef __AVX2__
inline void transpose8x8_ps(float *out, float *in, int m, int n)
{
	__m256 row0 = _mm256_load_ps(in+0*n);
	__m256 row1 = _mm256_load_ps(in+1*n);
	__m256 row2 = _mm256_load_ps(in+2*n);
	__m256 row3 = _mm256_load_ps(in+3*n);
	__m256 row4 = _mm256_load_ps(in+4*n);
	__m256 row5 = _mm256_load_ps(in+5*n);
	__m256 row6 = _mm256_load_ps(in+6*n);
	__m256 row7 = _mm256_load_ps(in+7*n);

	__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
	__m256 __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
	__t0 = _mm256_unpacklo_ps(row0, row1);
	__t1 = _mm256_unpackhi_ps(row0, row1);
	__t2 = _mm256_unpacklo_ps(row2, row3);
	__t3 = _mm256_unpackhi_ps(row2, row3);
	__t4 = _mm256_unpacklo_ps(row4, row5);
	__t5 = _mm256_unpackhi_ps(row4, row5);
	__t6 = _mm256_unpacklo_ps(row6, row7);
	__t7 = _mm256_unpackhi_ps(row6, row7);
	__tt0 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
	__tt1 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
	__tt2 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
	__tt3 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
	__tt4 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
	__tt5 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
	__tt6 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
	__tt7 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
	row0 = _mm256_permute2f128_ps(__tt0, __tt4, 0x20);
	row1 = _mm256_permute2f128_ps(__tt1, __tt5, 0x20);
	row2 = _mm256_permute2f128_ps(__tt2, __tt6, 0x20);
	row3 = _mm256_permute2f128_ps(__tt3, __tt7, 0x20);
	row4 = _mm256_permute2f128_ps(__tt0, __tt4, 0x31);
	row5 = _mm256_permute2f128_ps(__tt1, __tt5, 0x31);
	row6 = _mm256_permute2f128_ps(__tt2, __tt6, 0x31);
	row7 = _mm256_permute2f128_ps(__tt3, __tt7, 0x31);

	_mm256_store_ps(out+0*m, row0);
	_mm256_store_ps(out+1*m, row1);
	_mm256_store_ps(out+2*m, row2);
	_mm256_store_ps(out+3*m, row3);
	_mm256_store_ps(out+4*m, row4);
	_mm256_store_ps(out+5*m, row5);
	_mm256_store_ps(out+6*m, row6);
	_mm256_store_ps(out+7*m, row7);
}
void transpose_AVX2(float *out, float *in, int m, int n)
{
	assert(m%8 == 0);
	assert(n%8 == 0);

	const int tilex = 64;
	const int tiley = 64;

	int blockx = n/tilex;
	int blocky = m/tiley;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (int i=0; i<m; i+=tiley)
	{
		for (int j=0; j<n; j+=tilex)
		{
			int mm = std::min(i+tiley, m);
			int nn = std::min(j+tilex, n);

			for (int ii=i; ii<mm; ii+=8)
			{
				for (int jj=j; jj<nn; jj+=8)
				{
					transpose8x8_ps(out+(jj*m+ii), in+(ii*n+jj), m, n);
				}
			}
		}
	}
}
#endif

void runMedian(float *data, float *datMedian, long int size, int w)
{
	AVLTree<float> treeMedian;

	w = w > size ? size : w;
	long int k = 0;

	int a = -floor(w * 0.5);
	int b = ceil(w * 0.5);

	for (long int i = 0; i < b; i++)
	{
		treeMedian.insertValue(data[i]);
	}
	datMedian[k++] = treeMedian.getMedian();

	for (long int i = 1; i < size; i++)
	{
		a++;
		b++;
		if (a > 0)
			treeMedian.removeValue(data[a - 1]);
		if (b <= size)
			treeMedian.insertValue(data[b - 1]);

		datMedian[k++] = treeMedian.getMedian();
	}
}

template <typename T>
void runMedian2(T *data, T *datMedian, long int size, int w)
{
	w = w > size ? size : w;

	std::multiset<T, greater<T>> lowhalf;
	std::multiset<T, less<T>> highhalf;

	int a = 0-w/2-1;
	int b = 0+(w-1)/2;

	lowhalf.insert(data[0]);
	T median = data[0];

	for (long int i=1; i<b; i++)
	{
		if (data[i] >= median)
		{
			highhalf.insert(data[i]);
		}
		else
		{
			lowhalf.insert(data[i]);
		}

		if (lowhalf.size() > highhalf.size()+1)
		{
			highhalf.insert(*lowhalf.begin());
			lowhalf.erase(lowhalf.begin());
		}
		else if (highhalf.size() > lowhalf.size()+1)
		{
			lowhalf.insert(*highhalf.begin());
			highhalf.erase(highhalf.begin());
		}

		if (lowhalf.size() > highhalf.size())
		{
			median = *lowhalf.begin();
		}
		else if (highhalf.size() > lowhalf.size())
		{
			median = *highhalf.begin();
		}
		else
		{
			median = (*lowhalf.begin()+*highhalf.begin())/2;
		}
	}

	for (long int i=0; i<w/2+1; i++)
	{
		if (data[b] >= median)
		{
			highhalf.insert(data[b]);
		}
		else
		{
			lowhalf.insert(data[b]);
		}

		if (lowhalf.size() > highhalf.size()+1)
		{
			highhalf.insert(*lowhalf.begin());
			lowhalf.erase(lowhalf.begin());
		}
		else if (highhalf.size() > lowhalf.size()+1)
		{
			lowhalf.insert(*highhalf.begin());
			highhalf.erase(highhalf.begin());
		}

		if (lowhalf.size() > highhalf.size())
		{
			median = *lowhalf.begin();
		}
		else if (highhalf.size() > lowhalf.size())
		{
			median = *highhalf.begin();
		}
		else
		{
			median = (*lowhalf.begin()+*highhalf.begin())/2;
		}

		datMedian[i] = median;

		a++;
		b++;
	}

	for (long int i=w/2+1; i<size-(w-1)/2; i++)
	{
		if (data[b] >= median)
		{
			highhalf.insert(data[b]);
		}
		else
		{
			lowhalf.insert(data[b]);
		}

		auto it = lowhalf.find(data[a]);
		if (it != lowhalf.end())
			lowhalf.erase(it);
		else
			highhalf.erase(highhalf.find(data[a]));

		if (lowhalf.size() > highhalf.size()+1)
		{
			highhalf.insert(*lowhalf.begin());
			lowhalf.erase(lowhalf.begin());
		}
		else if (highhalf.size() > lowhalf.size()+1)
		{
			lowhalf.insert(*highhalf.begin());
			highhalf.erase(highhalf.begin());
		}

		if (lowhalf.size() > highhalf.size())
		{
			median = *lowhalf.begin();
		}
		else if (highhalf.size() > lowhalf.size())
		{
			median = *highhalf.begin();
		}
		else
		{
			median = (*lowhalf.begin()+*highhalf.begin())/2;
		}

		datMedian[i] = median;

		a++;
		b++;
	}

	for (long int i=size-(w-1)/2; i<size; i++)
	{
		auto it = lowhalf.find(data[a]);
		if (it != lowhalf.end())
			lowhalf.erase(it);
		else
			highhalf.erase(highhalf.find(data[a]));

		if (lowhalf.size() > highhalf.size()+1)
		{
			highhalf.insert(*lowhalf.begin());
			lowhalf.erase(lowhalf.begin());
		}
		else if (highhalf.size() > lowhalf.size()+1)
		{
			lowhalf.insert(*highhalf.begin());
			highhalf.erase(highhalf.begin());
		}

		if (lowhalf.size() > highhalf.size())
		{
			median = *lowhalf.begin();
		}
		else if (highhalf.size() > lowhalf.size())
		{
			median = *highhalf.begin();
		}
		else
		{
			median = (*lowhalf.begin()+*highhalf.begin())/2;
		}

		datMedian[i] = median;

		a++;
		b++;
	}
}

template <typename T>
void runMedian3(T *data, T *datMedian, long int size, int w)
{
	multiset<T, greater<T>> treeMedian;

	w = w > size ? size : w;
	long int k = 0;

	int a = -floor(w * 0.5);
	int b = ceil(w * 0.5);

	for (long int i = 0; i < b; i++)
	{
		treeMedian.insert(data[i]);
	}

	if (treeMedian.size()%2 == 1)
	{
		int tmp = 0;
		for (auto it=treeMedian.begin(); it!=treeMedian.end(); ++it)
		{
			if (tmp++ == treeMedian.size()/2)
			{
				datMedian[k++] = *it;
			}
		}
	}
	else
	{
		int tmp = 0;
		for (auto it=treeMedian.begin(); it!=treeMedian.end(); ++it)
		{
			if (tmp == treeMedian.size()/2-1)
			{
				datMedian[k] = *it;
			}

			if (tmp == treeMedian.size()/2)
			{
				datMedian[k] += *it;
			}
			tmp++;
		}
		datMedian[k] /= 2;
		k++;
	}

	for (long int i = 1; i < size; i++)
	{
		a++;
		b++;
		if (a > 0)
			treeMedian.erase(treeMedian.find(data[a - 1]));
		if (b <= size)
			treeMedian.insert(data[b - 1]);

		if (treeMedian.size()%2 == 1)
		{
			int tmp = 0;
			for (auto it=treeMedian.begin(); it!=treeMedian.end(); ++it)
			{
				if (tmp++ == treeMedian.size()/2)
				{
					datMedian[k++] = *it;
				}
			}
		}
		else
		{
			int tmp = 0;
			for (auto it=treeMedian.begin(); it!=treeMedian.end(); ++it)
			{
				if (tmp == treeMedian.size()/2-1)
				{
					datMedian[k] = *it;
				}

				if (tmp == treeMedian.size()/2)
				{
					datMedian[k] += *it;
				}
				tmp++;
			}
			datMedian[k] /= 2;
			k++;
		}
	}
}

void cmul(vector<complex<float>> &x, vector<complex<float>> &y)
{
	assert(x.size() == y.size());

	long int size = x.size();

	vector<float> xr(size);
	vector<float> xi(size);
	vector<float> yr(size);
	vector<float> yi(size);
	
	for (long int i=0; i<size; i++)
	{
		xr[i] = real(x[i]);
		xi[i] = imag(x[i]);
	}
	for (long int i=0; i<size; i++)
	{
		yr[i] = real(y[i]);
		yi[i] = imag(y[i]);
	}

	for (long int i=0; i<size; i++)
	{
		float real = xr[i]*yr[i]-xi[i]*yi[i];
		float imag = xr[i]*yi[i]+xi[i]*yr[i];
		x[i].real(real);
		x[i].imag(imag);
	}
}

void get_s_radec(double ra, double dec, string &s_ra, string &s_dec)
{
	int ra_hh, ra_mm;
	float ra_ss;

	int ra_sign = ra>=0 ? 1:-1;
	ra *= ra_sign;
	ra_hh = ra/10000;
	ra_mm = (ra-ra_hh*10000)/100;
	ra_ss = ra-ra_hh*10000-ra_mm*100;

	int dec_dd, dec_mm;
	float dec_ss;

	int dec_sign = dec>=0 ? 1:-1;
	dec *= dec_sign;
	dec_dd = dec/10000;
	dec_mm = (dec-dec_dd*10000)/100;
	dec_ss = dec-dec_dd*10000-dec_mm*100;

	stringstream ss_ra_hh;
	ss_ra_hh << setw(2) << setfill('0') << ra_hh;
	string s_ra_hh = ss_ra_hh.str();

	stringstream ss_ra_mm;
	ss_ra_mm << setw(2) << setfill('0') << ra_mm;
	string s_ra_mm = ss_ra_mm.str();

	stringstream ss_ra_ss;
	ss_ra_ss << setprecision(2) << setw(5) << setfill('0') << fixed << ra_ss;
	string s_ra_ss = ss_ra_ss.str();

	stringstream ss_dec_dd;
	ss_dec_dd << setw(2) << setfill('0') << dec_dd;
	string s_dec_dd = ss_dec_dd.str();

	stringstream ss_dec_mm;
	ss_dec_mm << setw(2) << setfill('0') << dec_mm;
	string s_dec_mm = ss_dec_mm.str();

	stringstream ss_dec_ss;
	ss_dec_ss << setprecision(2) << setw(5) << setfill('0') << fixed << dec_ss;
	string s_dec_ss = ss_dec_ss.str();

	s_ra = s_ra_hh + ":" + s_ra_mm + ":" + s_ra_ss;
	if (dec_sign < 0)
		s_dec = "-" + s_dec_dd + ":" + s_dec_mm + ":" + s_dec_ss;
	else
		s_dec = s_dec_dd + ":" + s_dec_mm + ":" + s_dec_ss;
}

bool inverse_matrix4x4(const double m[16], double invOut[16])
{
	double inv[16], det;
	int i;

	inv[0] = m[5]  * m[10] * m[15] - 
			 m[5]  * m[11] * m[14] - 
			 m[9]  * m[6]  * m[15] + 
			 m[9]  * m[7]  * m[14] +
			 m[13] * m[6]  * m[11] - 
			 m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] + 
			  m[4]  * m[11] * m[14] + 
			  m[8]  * m[6]  * m[15] - 
			  m[8]  * m[7]  * m[14] - 
			  m[12] * m[6]  * m[11] + 
			  m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] - 
			 m[4]  * m[11] * m[13] - 
			 m[8]  * m[5] * m[15] + 
			 m[8]  * m[7] * m[13] + 
			 m[12] * m[5] * m[11] - 
			 m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] + 
			   m[4]  * m[10] * m[13] +
			   m[8]  * m[5] * m[14] - 
			   m[8]  * m[6] * m[13] - 
			   m[12] * m[5] * m[10] + 
			   m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] + 
			  m[1]  * m[11] * m[14] + 
			  m[9]  * m[2] * m[15] - 
			  m[9]  * m[3] * m[14] - 
			  m[13] * m[2] * m[11] + 
			  m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] - 
			 m[0]  * m[11] * m[14] - 
			 m[8]  * m[2] * m[15] + 
			 m[8]  * m[3] * m[14] + 
			 m[12] * m[2] * m[11] - 
			 m[12] * m[3] * m[10];

	inv[9] = -m[0]  * m[9] * m[15] + 
			  m[0]  * m[11] * m[13] + 
			  m[8]  * m[1] * m[15] - 
			  m[8]  * m[3] * m[13] - 
			  m[12] * m[1] * m[11] + 
			  m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] - 
			  m[0]  * m[10] * m[13] - 
			  m[8]  * m[1] * m[14] + 
			  m[8]  * m[2] * m[13] + 
			  m[12] * m[1] * m[10] - 
			  m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] - 
			 m[1]  * m[7] * m[14] - 
			 m[5]  * m[2] * m[15] + 
			 m[5]  * m[3] * m[14] + 
			 m[13] * m[2] * m[7] - 
			 m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] + 
			  m[0]  * m[7] * m[14] + 
			  m[4]  * m[2] * m[15] - 
			  m[4]  * m[3] * m[14] - 
			  m[12] * m[2] * m[7] + 
			  m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] - 
			  m[0]  * m[7] * m[13] - 
			  m[4]  * m[1] * m[15] + 
			  m[4]  * m[3] * m[13] + 
			  m[12] * m[1] * m[7] - 
			  m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] + 
			   m[0]  * m[6] * m[13] + 
			   m[4]  * m[1] * m[14] - 
			   m[4]  * m[2] * m[13] - 
			   m[12] * m[1] * m[6] + 
			   m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] + 
			  m[1] * m[7] * m[10] + 
			  m[5] * m[2] * m[11] - 
			  m[5] * m[3] * m[10] - 
			  m[9] * m[2] * m[7] + 
			  m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] - 
			 m[0] * m[7] * m[10] - 
			 m[4] * m[2] * m[11] + 
			 m[4] * m[3] * m[10] + 
			 m[8] * m[2] * m[7] - 
			 m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] + 
			   m[0] * m[7] * m[9] + 
			   m[4] * m[1] * m[11] - 
			   m[4] * m[3] * m[9] - 
			   m[8] * m[1] * m[7] + 
			   m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] - 
			  m[0] * m[6] * m[9] - 
			  m[4] * m[1] * m[10] + 
			  m[4] * m[2] * m[9] + 
			  m[8] * m[1] * m[6] - 
			  m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}

bool get_inverse_matrix3x3(const double m[3][3], double invOut[3][3])
{
	double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
				m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
				m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

	if (det == 0) return false;

	double invdet = 1 / det;

	invOut[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
	invOut[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
	invOut[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
	invOut[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
	invOut[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
	invOut[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
	invOut[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
	invOut[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
	invOut[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

	return true;
}

template <typename T>
bool get_error_from_chisq_matrix(T &xerr, T &yerr, vector<T> &x, vector<T> &y, vector<T> &mxchisq)
{
	long int m = y.size();
	long int n = x.size();
	assert(mxchisq.size() == m*n);

	double A[4][4];
	double B[4];
	for (long int i=0; i<4; i++)
	{   
		B[i] = 0.;
		for (long int j=0; j<4; j++)
		{
			A[i][j] = 0.;
		}
	}

	vector<T> x2(n), x3(n), x4(n);
	vector<T> y2(m), y3(m), y4(m);

	for (long int j=0; j<n; j++)
	{
		x2[j] = x[j]*x[j];
		x3[j] = x2[j]*x[j];
		x4[j] = x3[j]*x[j];
	}

	for (long int i=0; i<m; i++)
	{
		y2[i] = y[i]*y[i];
		y3[i] = y2[i]*y[i];
		y4[i] = y3[i]*y[i];
	}    

	for (long int i=0; i<m; i++)
	{
		for (long int j=0; j<n; j++)
		{
			A[0][0] += x4[j];
			A[0][1] = A[1][0] += 2.*x3[j]*y[i];
			A[0][2] = A[2][0] += x2[j]*y2[i];
			A[0][3] = A[3][0] += x2[j];
			A[1][1] += 4.*x2[j]*y2[i];
			A[1][2] = A[2][1] += 2*x[j]*y3[i];
			A[1][3] = A[3][1] += 2*x[j]*y[i];
			A[2][2] += y4[i];
			A[2][3] = A[3][2] += y2[i];
			A[3][3] += 1.;
			B[0] += x2[j]*mxchisq[i*n+j];
			B[1] += 2.*x[j]*y[i]*mxchisq[i*n+j];
			B[2] += y2[i]*mxchisq[i*n+j];
			B[3] += mxchisq[i*n+j];
		}
	}

	double invA[4][4];

	if (!get_inverse_matrix4x4(A, invA)) return false;

	double a=0., b=0., c=0., d=0.;

	for (long int i=0; i<4; i++)
	{
		a += invA[0][i]*B[i];
		b += invA[1][i]*B[i];
		c += invA[2][i]*B[i];
		d += invA[3][i]*B[i];
	}

	xerr = sqrt(abs(c/(a*c-b*b)));
	yerr = sqrt(abs(a/(a*c-b*b)));

	return true;
}

template <typename T>
bool get_error_from_chisq_matrix(T &xerr, vector<T> &x, vector<T> &vchisq)
{
	int n = x.size();
	assert(n == vchisq.size());

	double A[3][3];
	double B[3];
	for (long int i=0; i<3; i++)
	{
		B[i] = 0.;
		for (long int j=0; j<3; j++)
		{
			A[i][j] = 0.;
		}
	}

	vector<T> x2(n), x3(n), x4(n);
	for (long int i=0; i<n; i++)
	{
		x2[i] = x[i]*x[i];
		x3[i] = x2[i]*x[i];
		x4[i] = x3[i]*x[i];
	}

	for (long int i=0; i<n; i++)
	{
		A[0][0] += x4[i];
		A[0][1] = A[1][0] += x3[i];
		A[0][2] = A[2][0] += x2[i];
		A[1][1] += x2[i];
		A[1][2] = A[2][1] += x[i];
		A[2][2] += 1.;
		
		B[0] += x2[i]*vchisq[i];
		B[1] += x[i]*vchisq[i];
		B[2] += vchisq[i];
	}

	double invA[3][3];

	if (!get_inverse_matrix3x3(A, invA)) return false;

	double a=0., b=0., c=0.;
	for (long int i=0; i<3; i++)
	{
		a += invA[0][i]*B[i];
		b += invA[1][i]*B[i];
		c += invA[2][i]*B[i];
	}

	xerr = sqrt(abs(1./a));

	return true;
}

#ifdef HAVE_SOFA
#include "sofa.h"
/**
 * @brief calculate source GB and GL using sofa lib 
 * 
 * @param gl : deg
 * @param gb : deg
 * @param s_ra: 00:00:00
 * @param s_dec: 00:00:00
 */
void get_gl_gb(double &gl, double &gb, const std::string &s_ra, const std::string &s_dec)
{
	double ra=0., dec=0.;
	get_rad_radec(s_ra, s_dec, ra, dec);
	iauIcrs2g(ra, dec, &gl, &gb);
	gl *= 180./M_PI;
	gb *= 180./M_PI;
}
#endif

#ifdef HAVE_YMW16
#include "cn.h"
/**
 * @brief Get the maxdm ymw16 lib
 * 
 * @param gl: degree
 * @param gb: degree
 * @return maxdm : pc/cc
 */
double get_maxdm_ymw16(double gl, double gb)
{
	double ymw16_maxdm = 0.;
	
	if (std::getenv("YMW16_DIR") == NULL)
	{
		std::cerr<<"Warning: environment variable YMW16_DIR not set. DM YMW16 and Distance YMW16 will not be calculated."<<endl;
	}
	else
	{
		char dirname[1024];
		std::strcpy(dirname, std::getenv("YMW16_DIR"));
		char text[1024]="\0";

		ymw16_maxdm = dmdtau(gl, gb, 1e6, 0, 2, 1, 0, dirname, text);    
	}
	
	return ymw16_maxdm;
}
/**
 * @brief Get the dist ymw16 lib
 * 
 * @param gl : degree
 * @param gb : degree
 * @param dm : pc/cc
 * @return dist : pc
 */
double get_dist_ymw16(double gl, double gb, double dm)
{
	double ymw16_dist = 0.;
	
	if (std::getenv("YMW16_DIR") == NULL)
	{
		std::cerr<<"Warning: environment variable YMW16_DIR not set. DM YMW16 and Distance YMW16 will not be calculated."<<endl;
	}
	else
	{
		char dirname[1024];
		std::strcpy(dirname, std::getenv("YMW16_DIR"));
		char text[1024]="\0";

		ymw16_dist = dmdtau(gl, gb, dm, 0, 1, 1, 0, dirname, text);
	}
	
	return ymw16_dist;
}
#endif

template <typename T>
void get_skewness_kurtosis(T profile, int size, double &skewness, double &kurtosis)
{
	long int nbin = size;

	double boxsum = 0.;
	for (long int i=0; i<nbin/2; i++)
	{
		boxsum += profile[i];
	}
	double min = boxsum;
	long int istart = 0;
	long int iend = nbin/2;
	for (long int i=0; i<nbin; i++)
	{
		boxsum -= profile[i];
		boxsum += profile[(i+nbin/2)%nbin];
		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/2+i+1;
		}
	}

	double tmp_mean1 = 0.;
	double tmp_mean2 = 0.;
	double tmp_mean3 = 0.;
	double tmp_mean4 = 0.;
	for (long int i=istart; i<iend; i++)
	{
		double tmp1 = profile[i%nbin];
		double tmp2 = tmp1*tmp1;
		double tmp3 = tmp2*tmp1;
		double tmp4 = tmp2*tmp2;
		tmp_mean1 += tmp1;
		tmp_mean2 += tmp2;
		tmp_mean3 += tmp3;
		tmp_mean4 += tmp4;
	}
	tmp_mean1 /= (nbin/2);
	tmp_mean2 /= (nbin/2);
	tmp_mean3 /= (nbin/2);
	tmp_mean4 /= (nbin/2);

	double tmp = tmp_mean1*tmp_mean1;
	double tmp_std = tmp_mean2-tmp;
	if (tmp_std == 0.)
	{
		skewness = 0.;
		kurtosis = 0.;
	}
	else
	{
		skewness = (tmp_mean3-3.*tmp_mean2*tmp_mean1+2.*tmp*tmp_mean1)/(tmp_std*std::sqrt(tmp_std));
		kurtosis = (tmp_mean4-4.*tmp_mean3*tmp_mean1+6.*tmp_mean2*tmp-3.*tmp*tmp)/(tmp_std*tmp_std) - 3.;
	}
}

template <typename T>
void get_mean_var(T profile, int size, double &mean, double &var)
{
	long int nbin = size;

	double boxsum = 0.;
	for (long int i=0; i<nbin/2; i++)
	{
		boxsum += profile[i];
	}
	double min = boxsum;
	long int istart = 0;
	long int iend = nbin/2;
	for (long int i=0; i<nbin; i++)
	{
		boxsum -= profile[i];
		boxsum += profile[(i+nbin/2)%nbin];
		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/2+i+1;
		}
	}

	double tmp_mean = min/(nbin/2);
	double tmp_var = 0.;
	for (long int i=istart; i<iend; i++)
	{
		double tmp = profile[i%nbin];
		tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
	}
	tmp_var /= (nbin/2);

	mean = tmp_mean;
	var = tmp_var;
}

template <typename T>
void get_mean_var2(T profile, int size, double &mean, double &var)
{
	long int nbin = size;

	double boxsum = 0.;
	for (long int i=0; i<nbin/4; i++)
	{
		boxsum += profile[i];
	}
	double min = boxsum;
	long int istart = 0;
	long int iend = nbin/4;
	for (long int i=0; i<nbin-1; i++)
	{
		boxsum -= profile[i];
		boxsum += profile[(i+nbin/4)%nbin];
		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/4+i+1;
		}
	}

	long int istart1 = istart;
	long int iend1 = iend;

	boxsum = 0.;
	for (long int i=iend1; i<iend1+nbin/4; i++)
	{
		boxsum += profile[i%nbin];
	}
	min = boxsum;
	istart = iend;
	iend = istart+nbin/4;
	for (long int i=iend1; i<iend1+nbin/2; i++)
	{
		boxsum -= profile[i%nbin];
		boxsum += profile[(i+nbin/4)%nbin];
		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/4+i+1;
		}
	}
	for (long int i=iend1+nbin/2; i<iend1+3*nbin/4-1; i++)
	{
		boxsum -= profile[i%nbin];
		boxsum += profile[(i+nbin/2)%nbin];
		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/2+i+1;
		}
	}

	long int istart2 = istart;
	long int iend2 = iend;

	double tmp_mean = 0.;
	double tmp_var = 0.;
	if (iend2-istart2 == nbin/2)
	{
		for (long int i=istart2; i<iend2; i++)
		{
			double tmp = profile[i%nbin];
			tmp_mean += tmp;
			tmp_var += tmp*tmp;
		}
	}
	else
	{
		for (long int i=istart1; i<iend1; i++)
		{
			double tmp = profile[i%nbin];
			tmp_mean += tmp;
			tmp_var += tmp*tmp;
		}
		for (long int i=istart2; i<iend2; i++)
		{
			double tmp = profile[i%nbin];
			tmp_mean += tmp;
			tmp_var += tmp*tmp;
		}
	}

	tmp_mean /= (nbin/2);
	tmp_var /= (nbin/2);
	tmp_var -= tmp_mean*tmp_mean;

	mean = tmp_mean;
	var = tmp_var;
}

template <typename T>
void get_mean_var(T profiles, int nrow, int ncol, double &mean, double &var)
{
	long int nchan = nrow;
	long int nbin = ncol;

	double boxsum = 0.;
	for (long int j=0; j<nchan; j++)
	{
		for (long int i=0; i<nbin/2; i++)
		{
			boxsum += profiles[j*nbin+i];
		}
	}
	double min = boxsum;
	long int istart = 0;
	long int iend = nbin/2;

	for (long int i=0; i<nbin; i++)
	{
		for (long int j=0; j<nchan; j++)
		{
			boxsum -= profiles[j*nbin+i];
			boxsum += profiles[j*nbin+(i+nbin/2)%nbin];
		}

		if (boxsum < min)
		{
			min = boxsum;
			istart = i+1;
			iend = nbin/2+i+1;
		}
	}

	double tmp_mean = min/((nbin/2)*nchan);
	double tmp_var = 0.;
	for (long int j=0; j<nchan; j++)
	{
		for (long int i=istart; i<iend; i++)
		{
			double tmp = profiles[j*nbin+i%nbin];
			tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
		}
	}
	tmp_var /= (nbin/2)*nchan;

	mean = tmp_mean;
	var = tmp_var;
}

template <typename T>
void get_mean_var(T profile, T profiles, int nsubint, int nchan, int nbin, double &mean, double &var)
{
	std::vector<double> alpha(nsubint*nchan, 0.);
	std::vector<double> beta(nsubint*nchan, 0.);

	std::vector<double> profile_sort(nbin, 0.);

	double se = 0., ss = 0.;
	for (long int i=0; i<nbin; i++)
	{
		se += profile[i];
		ss += profile[i]*profile[i];

		profile_sort[i] = profile[i];
	}

	double temp = se*se-ss*nbin;
	if (temp == 0) temp = 1.;

	for (long int k=0; k<nsubint; k++)
	{
		for (long int j=0; j<nchan; j++)
		{
			double xe = 0., xs = 0.;
			for (long int i=0; i<nbin; i++)
			{
				xe += profiles[k*nchan*nbin+j*nbin+i];
				xs += profiles[k*nchan*nbin+j*nbin+i]*profile[i];
			}

			alpha[k*nchan+j] = (se*xe-xs*nbin)/temp;
			beta[k*nchan+j] = (xs*se-xe*ss)/temp;
		}
	}

	mean = 0.;
	var = 0.;
	for (long int k=0; k<nsubint; k++)
	{
		for (long int j=0; j<nchan; j++)
		{
			double tmp_mean = 0.;
			double tmp_var = 0.;
			for (long int i=0; i<nbin; i++)
			{
				double tmp = profiles[k*nchan*nbin+j*nbin+i]-alpha[k*nchan+j]*profile[i]-beta[k*nchan+j];
				tmp_mean += tmp;
				tmp_var += tmp*tmp;
			}
			tmp_mean /= nbin;
			tmp_var /= nbin;
			tmp_var -= tmp_mean*tmp_mean;
			var += tmp_var;
		}
	}
}

template bool get_error_from_chisq_matrix<float>(float &xerr, float &yerr, vector<float> &x, vector<float> &y, vector<float> &mxchisq);
template bool get_error_from_chisq_matrix<float>(float &xerr, vector<float> &x, vector<float> &vchisq);
template bool get_error_from_chisq_matrix<double>(double &xerr, double &yerr, vector<double> &x, vector<double> &y, vector<double> &mxchisq);
template bool get_error_from_chisq_matrix<double>(double &xerr, vector<double> &x, vector<double> &vchisq);

template void get_skewness_kurtosis<std::vector<float>::iterator>(std::vector<float>::iterator profile, int size, double &skewness, double &kurtosis);
template void get_skewness_kurtosis<std::vector<double>::iterator>(std::vector<double>::iterator profile, int size, double &skewness, double &kurtosis);

template void get_mean_var<std::vector<float>::iterator>(std::vector<float>::iterator profile, int size, double &mean, double &var);
template void get_mean_var2<std::vector<float>::iterator>(std::vector<float>::iterator profile, int size, double &mean, double &var);
template void get_mean_var<std::vector<float>::iterator>(std::vector<float>::iterator profiles, int nrow, int ncol, double &mean, double &var);
template void get_mean_var<std::vector<double>::iterator>(std::vector<double>::iterator profile, int size, double &mean, double &var);
template void get_mean_var2<std::vector<double>::iterator>(std::vector<double>::iterator profile, int size, double &mean, double &var);
template void get_mean_var<std::vector<double>::iterator>(std::vector<double>::iterator profiles, int nrow, int ncol, double &mean, double &var);

template void get_mean_var<std::vector<double>::iterator>(std::vector<double>::iterator profile, std::vector<double>::iterator profiles, int nsubint, int nchan, int nbin, double &mean, double &var);
template void get_mean_var<std::vector<float>::iterator>(std::vector<float>::iterator profile, std::vector<float>::iterator profiles, int nsubint, int nchan, int nbin, double &mean, double &var);

template void transpose<float>(float *out, float *in, int m, int n);
template void transpose_pad<float>(float *out, float *in, int m, int n);
template void transpose_pad<float>(float *out, float *in, int m, int n, int tiley, int tilex);
template void transpose<complex<float>>(complex<float> *out, complex<float> *in, int m, int n);
template void transpose_pad<unsigned char>(unsigned char *out, unsigned char *in, int m, int n);
template void transpose_pad<unsigned char>(unsigned char *out, unsigned char *in, int m, int n, int tiley, int tilex);
template void transpose_pad<short>(short *out, short *in, int m, int n);
template void transpose_pad<short>(short *out, short *in, int m, int n, int tiley, int tilex);

template void transpose_pad<complex<float>>(complex<float> *out, complex<float> *in, int m, int n);
template void transpose_pad<complex<float>>(complex<float> *out, complex<float> *in, int m, int n, int tiley, int tilex);

template void transpose<double>(double *out, double *in, int m, int n);
template void transpose_pad<double>(double *out, double *in, int m, int n);
template void transpose_pad<double>(double *out, double *in, int m, int n, int tiley, int tilex);

template void runMedian2<float>(float *data, float *datMedian, long int size, int w);
template void runMedian2<double>(double *data, double *datMedian, long int size, int w);

template void runMedian3<float>(float *data, float *datMedian, long int size, int w);
