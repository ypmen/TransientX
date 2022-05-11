/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-23 20:50:35
 * @modify date 2021-04-23 20:50:35
 * @desc [description]
 */

#ifndef HOUGH_TRANSFORM_H
#define HOUGH_TRANSFORM_H


#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <assert.h>

void fadd2(float * x0, float * x1, float *y0, float *y1, int n)
{
	for (long int i=0; i<n; i++)
	{
		x1[i] = x0[i] + y1[i];
		x0[i] = x0[i] + y0[i];
	}
}

void transform( std::vector<float> &temp, std::vector<float> &cache0, std::vector<float> &cache1,  std::vector<float> &in, std::vector<std::vector<int>> &shift, int nrow, int ncol, int depth, int irow)
{
	int n = std::log2(nrow);

	if (depth == n)
	{
		std::memcpy(temp.data()+irow*ncol, in.data()+irow*ncol, sizeof(float)*ncol);

		return;
	}

	int nrowsub = nrow/pow(2, depth+1);

	transform(temp, cache0, cache1, in, shift, nrow, ncol, depth+1, 2*irow);
	transform(temp, cache0, cache1, in, shift, nrow, ncol, depth+1, 2*irow+1);

	for (long int k=0; k<nrowsub; k++)
	{
		int delay0 = shift[n-1-depth][irow*(nrowsub*2)+k];
		int delay1 = shift[n-1-depth][irow*(nrowsub*2)+(k+nrowsub)];
		
		std::rotate_copy(temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol, temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol+delay0, temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol+ncol, cache0.begin());
		std::rotate_copy(temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol, temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol+delay1, temp.begin()+(2*irow+1)*nrowsub*ncol+k*ncol+ncol, cache1.begin());

		fadd2(temp.data()+irow*(nrowsub*2)*ncol+k*ncol, temp.data()+irow*(nrowsub*2)*ncol+(k+nrowsub)*ncol, cache0.data(), cache1.data(), ncol);
	}
}

void get_shiftplan(std::vector<std::vector<int>> &shiftplan, std::vector<int> &hrange, int nrow, int ncol, std::function<int(int, int, int)> &get_shift)
{
	/* calculate shift */
	int nrow_tmp = nrow;
	int nrowsub_tmp = nrow/nrow_tmp;

	hrange.resize(nrow, 0);
	for (long int k=0; k<nrowsub_tmp; k++)
	{
		hrange[k] = k*nrow_tmp;
	}

	std::vector<int> xrange_tmp(nrow, 0.);
	for (long int j=0; j<nrow; j++)
	{
		xrange_tmp[j] = j;
	}

	int n = log2(nrow);
	for (long int l=0; l<n; l++)
	{
		std::vector<int> shift_tmp(nrow_tmp*nrowsub_tmp);
		for (long int j=0; j<nrow_tmp/2; j++)
		{
			for (long int k=0; k<nrowsub_tmp; k++)
			{
				hrange[k] = hrange[k];
				hrange[k+nrowsub_tmp] = hrange[k]+nrow_tmp/2;
			}

			for (long int k=0; k<nrowsub_tmp*2; k++)
			{
				int s = get_shift(xrange_tmp[j*2], xrange_tmp[j*2+1], hrange[k])%ncol;
				if (s < 0) s += ncol;
				
				shift_tmp[j*(nrowsub_tmp*2)+k] = s;
			}

			xrange_tmp[j] = xrange_tmp[j*2];
		}
		shiftplan.push_back(shift_tmp);

		nrowsub_tmp <<= 1;
		nrow_tmp >>= 1;
	}
}

void hough_transform(std::vector<float> &in, std::vector<float> &out, int nrow, int ncol, std::function<int(int, int, int)> &get_shift)
{
	std::vector<std::vector<int>> shiftplan;
	std::vector<int> hrange;
	get_shiftplan(shiftplan, hrange, nrow, ncol, get_shift);

	std::vector<float> temp(nrow*ncol, 0.);
	std::vector<float> cache0(ncol, 0.);
	std::vector<float> cache1(ncol, 0.);

	transform(temp, cache0, cache1, in, shiftplan, nrow, ncol, 0, 0);

	out.resize(nrow*ncol, 0.);
	for (long int j=0; j<nrow; j++)
	{
		int irow = hrange[j];
		for (long int i=0; i<ncol; i++)
		{
			out[irow*ncol+i] =  temp[j*ncol+i];
		}
	}
}


#endif /* HOUGH_TRANSFORM_H */
