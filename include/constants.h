/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-05-19 22:52:10
 * @modify date 2020-05-19 22:52:10
 * @desc [description]
 */

#ifndef CONSTANTS
#define CONSTANTS

#define CONST_C 299792458 /*m/s*/

inline double fdot2acc(double fdot, double f)
{
	return -fdot/f*CONST_C;
}

inline double acc2fdot(double acc, double f)
{
	return -f*acc/CONST_C;
}

#endif /* CONSTANTS */
