/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-05-14 10:30:24
 * @modify date 2022-05-14 10:30:24
 * @desc [description]
 */

#ifndef AVX2_H
#define AVX2_H

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <cmath>
#include <immintrin.h>

#ifndef __FMA__
inline __m256d mm256_fmadd_pd (__m256d __A, __m256d __B, __m256d __C)
{
	return _mm256_add_pd(_mm256_mul_pd(__A, __B), __C);
}

inline __m256 mm256_fmadd_ps (__m256 __A, __m256 __B, __m256 __C)
{
	return _mm256_add_ps(_mm256_mul_ps(__A, __B), __C);
}
#else
inline __m256d mm256_fmadd_pd (__m256d __A, __m256d __B, __m256d __C)
{
	return _mm256_fmadd_pd (__A, __B, __C);
}

inline __m256 mm256_fmadd_ps (__m256 __A, __m256 __B, __m256 __C)
{
	return _mm256_fmadd_ps ( __A, __B, __C);
}
#endif

namespace PulsarX {

typedef __attribute__(( aligned(32))) float aligned_float;
typedef __attribute__(( aligned(32))) double aligned_double;
typedef __attribute__(( aligned(32))) int aligned_int;
typedef __attribute__(( aligned(32))) unsigned char aligned_uchar;

inline void print(__m256 a)
{
	float *b = (float *)&a;
	std::cout<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<" "<<b[4]<<" "<<b[5]<<" "<<b[6]<<" "<<b[7]<<std::endl;
}

inline void print(__m256d a)
{
	double *b = (double *)&a;
	std::cout<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<std::endl;
}

inline void print(__m256i a)
{
	int *b = (int *)&a;
	std::cout<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<" "<<b[4]<<" "<<b[5]<<" "<<b[6]<<" "<<b[7]<<std::endl;
}

inline float hadd(__m256 a)
{
	__m256 b = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	a = _mm256_hadd_ps(b, a);
	a = _mm256_hadd_ps(b, a);

	return ((float *)&a)[3] + ((float *)&a)[7];
}

inline double haddd(__m256d a)
{
	__m256d b = _mm256_set_pd(0., 0., 0., 0.);
	a = _mm256_hadd_pd(a, b);

	return ((double *)&a)[0] + ((double *)&a)[2];
}

inline float reduce(const aligned_float * const data, size_t size)
{
	__m256 avx_acc = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	for (size_t i=0; i<size/8; i++)
	{
		__m256 axv_data = _mm256_load_ps(data+i*8);
		avx_acc = _mm256_add_ps(avx_acc, axv_data);
	}

	return hadd(avx_acc);
}

inline void accumulate_mean(
	aligned_double * const mean,
	aligned_double * const mean_scale,
	double scale,
	const aligned_float * const data,
	size_t size
)
{
	__m256d avx_scale = _mm256_set_pd(scale, scale, scale, scale);
	for (size_t i=0; i<size/4; i++)
	{
		__m128 avx_tmp = _mm_load_ps(data + i * 4);
		__m256d avx_data = _mm256_cvtps_pd(avx_tmp);
		__m256d avx_mean = _mm256_load_pd(mean + i * 4);
		__m256d avx_mean_scale = _mm256_load_pd(mean_scale + i * 4);

		avx_mean = _mm256_add_pd(avx_data, avx_mean);
		avx_mean_scale = mm256_fmadd_pd(avx_data, avx_scale, avx_mean_scale);

		_mm256_store_pd(mean + i * 4, avx_mean);
		_mm256_store_pd(mean_scale + i * 4, avx_mean_scale);
	}
}

inline void accumulate_mean_var(
	aligned_double * const mean,
	aligned_double * const var,
	const aligned_float * const data,
	size_t size
)
{
	for (size_t i=0; i<size/4; i++)
	{
		__m128 avx_tmp = _mm_load_ps(data + i * 4);
		__m256d avx_data = _mm256_cvtps_pd(avx_tmp);
		__m256d avx_mean = _mm256_load_pd(mean + i * 4);
		__m256d avx_var = _mm256_load_pd(var + i * 4);

		avx_mean = _mm256_add_pd(avx_mean, avx_data);
		avx_var = mm256_fmadd_pd(avx_data, avx_data, avx_var);

		_mm256_store_pd(mean + i * 4, avx_mean);
		_mm256_store_pd(var + i * 4, avx_var);
	}
}

inline void accumulate_mean_var_scale(
	aligned_double * const mean,
	aligned_double * const var,
	const aligned_float * const data,
	double scale,
	size_t size
)
{
	__m256d avx_scale = _mm256_set_pd(scale, scale, scale, scale);
	for (size_t i=0; i<size/4; i++)
	{
		__m128 avx_tmp = _mm_load_ps(data + i * 4);
		__m256d avx_data = _mm256_cvtps_pd(avx_tmp);
		avx_data = _mm256_mul_pd(avx_data, avx_scale);
		__m256d avx_mean = _mm256_load_pd(mean + i * 4);
		__m256d avx_var = _mm256_load_pd(var + i * 4);

		avx_mean = _mm256_add_pd(avx_mean, avx_data);
		avx_var = mm256_fmadd_pd(avx_data, avx_data, avx_var);

		_mm256_store_pd(mean + i * 4, avx_mean);
		_mm256_store_pd(var + i * 4, avx_var);
	}
}

inline void accumulate_mean_var2(
	double &mean,
	double &var,
	const aligned_float * const data,
	size_t size
)
{
	__m256d avx_mean = _mm256_set_pd(0., 0., 0., 0.);
	__m256d avx_var = _mm256_set_pd(0., 0., 0., 0.);
	for (size_t i=0; i<size/4; i++)
	{
		__m128 avx_tmp = _mm_load_ps(data + i * 4);
		__m256d avx_data = _mm256_cvtps_pd(avx_tmp);
		avx_mean = _mm256_add_pd(avx_data, avx_mean);
		avx_var = mm256_fmadd_pd(avx_data, avx_data, avx_var);
	}

	mean = haddd(avx_mean);
	var = haddd(avx_var);
}

inline void accumulate_mean_var3(
	float &mean,
	float &var,
	const aligned_float * const data,
	size_t size
)
{
	__m256 avx_mean = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	__m256 avx_var = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data = _mm256_load_ps(data + i * 8);
		avx_mean = _mm256_add_ps(avx_data, avx_mean);
		avx_var = mm256_fmadd_ps(avx_data, avx_data, avx_var);
	}

	mean = hadd(avx_mean);
	var = hadd(avx_var);
}

inline void accumulate_mean1_mean2_mean3_mean4_corr1(
	aligned_double * const mean1,
	aligned_double * const mean2,
	aligned_double * const mean3,
	aligned_double * const mean4,
	aligned_double * const corr1,
	aligned_double * const last_data,
	const aligned_float * const data,
	size_t size
)
{
	for (size_t i=0; i<size/4; i++)
	{
		__m128 avx_tmp = _mm_load_ps(data + i * 4);
		__m256d avx_data = _mm256_cvtps_pd(avx_tmp);

		__m256d avx_mean1 = _mm256_load_pd(mean1 + i * 4);
		__m256d avx_mean2 = _mm256_load_pd(mean2 + i * 4);
		__m256d avx_mean3 = _mm256_load_pd(mean3 + i * 4);
		__m256d avx_mean4 = _mm256_load_pd(mean4 + i * 4);
		__m256d avx_corr1 = _mm256_load_pd(corr1 + i * 4);
		__m256d avx_last_data = _mm256_load_pd(last_data + i * 4);

		__m256d avx_data2 = _mm256_mul_pd(avx_data, avx_data);

		avx_mean1 = _mm256_add_pd(avx_data, avx_mean1);
		avx_mean2 = mm256_fmadd_pd(avx_data, avx_data, avx_mean2);
		avx_mean3 = mm256_fmadd_pd(avx_data, avx_data2, avx_mean3);
		avx_mean4 = mm256_fmadd_pd(avx_data2, avx_data2, avx_mean4);
		avx_corr1 = mm256_fmadd_pd(avx_data, avx_last_data, avx_corr1);

		_mm256_store_pd(mean1 + i * 4, avx_mean1);
		_mm256_store_pd(mean2 + i * 4, avx_mean2);
		_mm256_store_pd(mean3 + i * 4, avx_mean3);
		_mm256_store_pd(mean4 + i * 4, avx_mean4);
		_mm256_store_pd(corr1 + i * 4, avx_corr1);
		_mm256_store_pd(last_data + i * 4, avx_data);
	}
}

inline void normalize(
	aligned_float * const data_out,
	const aligned_float * const data_in,
	const aligned_float * const weight,
	const aligned_float * const mean,
	const aligned_float * const stddev,
	size_t size
)
{
	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data_in = _mm256_load_ps(data_in + i * 8);
		__m256 avx_weight = _mm256_load_ps(weight + i * 8);
		__m256 avx_mean = _mm256_load_ps(mean + i * 8);
		__m256 avx_stddev = _mm256_load_ps(stddev + i * 8);

		__m256 avx_data_out = _mm256_mul_ps(_mm256_div_ps(_mm256_sub_ps(avx_data_in, avx_mean), avx_stddev), avx_weight);
		_mm256_store_ps(data_out + i * 8, avx_data_out);
	}
}

inline void normalize2(
	aligned_float * const data_out,
	const aligned_float * const data_in,
	const aligned_float * const mean,
	const aligned_float * const stddev_inv,
	size_t size
)
{
	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data_in = _mm256_load_ps(data_in + i * 8);

		__m256 avx_mean = _mm256_load_ps(mean + i * 8);
		__m256 avx_stddev_inv = _mm256_load_ps(stddev_inv + i * 8);

		__m256 avx_data_out = _mm256_mul_ps(_mm256_sub_ps(avx_data_in, avx_mean), avx_stddev_inv);
		_mm256_store_ps(data_out + i * 8, avx_data_out);
	}
}

inline void remove_baseline(
	aligned_float * const data_out,
	const aligned_float * const data_in,
	const aligned_float * const a,
	const aligned_float * const b,
	float s,
	size_t size
)
{
	__m256 avx_s = _mm256_set_ps(s, s, s, s, s, s, s, s);
	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data_in = _mm256_load_ps(data_in + i * 8);
		__m256 avx_a = _mm256_load_ps(a + i * 8);
		__m256 avx_b = _mm256_load_ps(b + i * 8);

		__m256 avx_data_out = _mm256_sub_ps(avx_data_in, _mm256_add_ps(_mm256_mul_ps(avx_a, avx_s), avx_b));

		_mm256_store_ps(data_out + i * 8, avx_data_out);
	}
}

inline float remove_baseline_reduce(
	aligned_float * const data_out,
	const aligned_float * const data_in,
	const aligned_float * const a,
	const aligned_float * const b,
	float s,
	size_t size
)
{
	__m256 avx_acc = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	__m256 avx_s = _mm256_set_ps(s, s, s, s, s, s, s, s);
	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data_in = _mm256_load_ps(data_in + i * 8);
		__m256 avx_a = _mm256_load_ps(a + i * 8);
		__m256 avx_b = _mm256_load_ps(b + i * 8);

		__m256 avx_data_out = _mm256_sub_ps(avx_data_in, _mm256_add_ps(_mm256_mul_ps(avx_a, avx_s), avx_b));

		_mm256_store_ps(data_out + i * 8, avx_data_out);

		avx_acc = mm256_fmadd_ps(avx_data_out, avx_data_out, avx_acc);
	}

	return hadd(avx_acc);
}

inline void kadane2D (
	aligned_float * const maxSum,
	aligned_int * const start,
	aligned_int * const end,
	const aligned_float * const arr,
	size_t nrow,
	size_t ncol
)
{
	float *sum = (float *) aligned_alloc(32, ncol*sizeof(float));
	int *local_start = (int *) aligned_alloc(32, ncol*sizeof(sizeof(int)));
	for (size_t j=0; j<ncol; j++)
	{
		sum[j] = 0.;
		maxSum[j] = -std::numeric_limits<float>::infinity();
		start[j] = 0;
		end[j] = -1;
		local_start[j] = 0;
	}

	__m256 avx_zero = _mm256_set_ps(0., 0., 0., 0., 0., 0., 0., 0.);
	__m256i avx_one = _mm256_set_epi32(1, 1, 1, 1, 1, 1, 1, 1);

	for (size_t i=0; i<nrow; i++)
	{
		__m256i avx_i = _mm256_set_epi32(i, i, i, i, i, i, i, i);
		__m256i avx_i1 = _mm256_add_epi32(avx_i, avx_one);

		for (size_t j=0; j<ncol/8; j++)
		{
			__m256 avx_arr = _mm256_load_ps(arr + i * ncol + j * 8);
			__m256 avx_sum = _mm256_load_ps(sum + j * 8);
			__m256 avx_maxSum = _mm256_load_ps(maxSum + j * 8);
			__m256i avx_local_start = _mm256_load_si256((__m256i *)(local_start + j * 8));
			__m256i avx_start = _mm256_load_si256((__m256i *)(start + j * 8));
			__m256i avx_end = _mm256_load_si256((__m256i *)(end + j * 8));

			avx_sum = _mm256_add_ps(avx_arr, avx_sum);

			__m256 avx_mask1 = _mm256_cmp_ps(avx_sum, avx_zero, 1);
			__m256 avx_mask2 = _mm256_cmp_ps(avx_sum, avx_zero, 13);
			__m256 avx_mask3 = _mm256_cmp_ps(avx_sum, avx_maxSum, 14);
			__m256 avx_mask4 = _mm256_and_ps(avx_mask3, avx_mask2);

			avx_sum = _mm256_blendv_ps(avx_sum, avx_zero, avx_mask1);
			avx_local_start = _mm256_blendv_epi8(avx_local_start, avx_i1, _mm256_castps_si256(avx_mask1));

			avx_maxSum = _mm256_blendv_ps(avx_maxSum, avx_sum, avx_mask4);
			avx_start = _mm256_blendv_epi8(avx_start, avx_local_start, _mm256_castps_si256(avx_mask4));
			avx_end = _mm256_blendv_epi8(avx_end, avx_i, _mm256_castps_si256(avx_mask4));
		
			_mm256_store_ps(sum + j * 8, avx_sum);
			_mm256_store_ps(maxSum + j * 8, avx_maxSum);
			_mm256_store_si256((__m256i *)(local_start + j * 8), avx_local_start);
			_mm256_store_si256((__m256i *)(start + j * 8), avx_start);
			_mm256_store_si256((__m256i *)(end + j * 8), avx_end);
		}
	}

	__m256i avx_zeroi = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
	__m256i avx_negone = _mm256_set_epi32(-1, -1, -1, -1, -1, -1, -1, -1);

	for (size_t j=0; j<ncol/8; j++)
	{
		__m256 avx_arr = _mm256_load_ps(arr + j * 8);
		__m256 avx_maxSum = _mm256_load_ps(maxSum + j * 8);
		__m256i avx_start = _mm256_load_si256((__m256i *)(start + j * 8));
		__m256i avx_end = _mm256_load_si256((__m256i *)(end + j * 8));
		
		__m256i avx_mask = _mm256_cmpeq_epi32(avx_end, avx_negone);
		avx_maxSum = _mm256_blendv_ps(avx_maxSum, avx_arr, _mm256_castsi256_ps(avx_mask));
		avx_start = _mm256_blendv_epi8(avx_start, avx_zeroi, avx_mask);
		avx_end = _mm256_blendv_epi8(avx_end, avx_zeroi, avx_mask);

		_mm256_store_ps(maxSum + j * 8, avx_maxSum);
		_mm256_store_si256((__m256i *)(start + j * 8), avx_start);
		_mm256_store_si256((__m256i *)(end + j * 8), avx_end);
	}

	free(sum);
	free(local_start);
}

inline void scale(
	aligned_uchar * const data_out,
	const aligned_float * const data_in,
	float scl,
	float offs,
	unsigned int nbits,
	size_t size
)
{
	__m256 avx_scl = _mm256_set_ps(scl, scl, scl, scl, scl, scl, scl, scl);
	__m256 avx_offs = _mm256_set_ps(offs, offs, offs, offs, offs, offs, offs, offs);
	float min = 0;
	float max = std::pow(2, nbits) - 1.;
	__m256 avx_min = _mm256_set_ps(min, min, min, min, min, min, min, min);
	__m256 avx_max = _mm256_set_ps(max, max, max, max, max, max, max, max);

	for (size_t i=0; i<size/8; i++)
	{
		__m256 avx_data = _mm256_load_ps(data_in + i * 8);
		avx_data = mm256_fmadd_ps(avx_data, avx_scl, avx_offs);
		avx_data = _mm256_round_ps(avx_data, 0);
		__m256 mask_min = _mm256_cmp_ps(avx_data, avx_min, 1);
		__m256 mask_max = _mm256_cmp_ps(avx_data, avx_max, 14);
		avx_data = _mm256_blendv_ps(avx_data, avx_min, mask_min);
		avx_data = _mm256_blendv_ps(avx_data, avx_max, mask_max);
		__m256i avx_data_out = _mm256_cvtps_epi32(avx_data);
		switch (nbits)
		{
		case 8:
		{
			int *tmp = (int *)&avx_data_out;
			data_out[i*8+0] = tmp[0];
			data_out[i*8+1] = tmp[1];
			data_out[i*8+2] = tmp[2];
			data_out[i*8+3] = tmp[3];
			data_out[i*8+4] = tmp[4];
			data_out[i*8+5] = tmp[5];
			data_out[i*8+6] = tmp[6];
			data_out[i*8+7] = tmp[7];
		};break;
		case 4:
		{
			__m256i avx_count = _mm256_set_epi32(4, 0, 4, 0, 4, 0, 4, 0);
			avx_data_out = _mm256_sllv_epi32(avx_data_out, avx_count);
			avx_data_out = _mm256_hadd_epi32(avx_data_out, avx_data_out);
			int *tmp = (int *)&avx_data_out;
			data_out[i*4+0] = tmp[0];
			data_out[i*4+1] = tmp[1];
			data_out[i*4+2] = tmp[4];
			data_out[i*4+3] = tmp[5];
		};break;
		case 2:
		{
			__m256i avx_count = _mm256_set_epi32(6, 4, 2, 0, 6, 4, 2, 0);
			avx_data_out = _mm256_sllv_epi32(avx_data_out, avx_count);
			avx_data_out = _mm256_hadd_epi32(avx_data_out, avx_data_out);
			avx_data_out = _mm256_hadd_epi32(avx_data_out, avx_data_out);
			int *tmp = (int *)&avx_data_out;
			data_out[i*2+0] = tmp[0];
			data_out[i*2+1] = tmp[4];
		};break;
		case 1:
		{
			__m256i avx_count = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
			avx_data_out = _mm256_sllv_epi32(avx_data_out, avx_count);
			avx_data_out = _mm256_hadd_epi32(avx_data_out, avx_data_out);
			avx_data_out = _mm256_hadd_epi32(avx_data_out, avx_data_out);
			int *tmp = (int *)&avx_data_out;
			data_out[i] = tmp[0] + tmp[4];
		};break;
		default:
		{
			std::cerr<<"Error: data type unsupported"<<std::endl;
			exit(-1);
		};break;
		}
	}
}

}


#endif /* AVX2_H */
