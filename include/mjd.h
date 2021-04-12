/*
 * mjd.h
 *
 *  Created on: Feb 25, 2020
 *      Author: ypmen
 */

#ifndef MJD_H_
#define MJD_H_

#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

class MJD
{
public:
	MJD(){stt_imjd = 0; stt_smjd = 0; stt_offs = 0.;}
	MJD(long int imjd, long int smjd, double offs){stt_imjd = imjd; stt_smjd = smjd; stt_offs = offs;}
	MJD(long double mjd){format(mjd);}
	~MJD(){}
	long double to_day(){return (long double)stt_imjd + (long double)stt_smjd/86400. + (long double)stt_offs/86400.;} const
	long double to_second(){return (long double)(stt_imjd*86400)+(long double)stt_smjd+(long double)stt_offs;} const
	bool to_date(char *dstr, int len, const char *format) const
	{
		struct tm greg;

		if (!gregorian (&greg, NULL))
		return false;

		if (stt_offs>=0.5)
			greg.tm_sec += 1;

		if (strftime (dstr, len, format, &greg) == 0)
		return false;

		return true;
	}
	void format()
	{
		long int nsec = floor(stt_offs);
		double fsec = stt_offs-nsec;
		stt_smjd += nsec;
		stt_offs = fsec;

		long int nday = stt_smjd/86400;
		nsec = stt_smjd%86400;
		if (nsec<0)
		{
			nsec += 86400;
			nday--;
		}

		stt_smjd = nsec;
		stt_imjd += nday;
	}
	void format(long double mjd)
	{
		stt_imjd = floor(mjd);
		long double sec = (mjd-(long double)stt_imjd)*(long double)86400.;
		stt_smjd = floor(sec);
		stt_offs = sec-stt_smjd;
	}
	MJD dividedby2()
	{
		long int nimjd = stt_imjd/2;
		long int nsec = (stt_imjd%2)*86400/2+stt_smjd/2;
		nimjd += nsec/86400;
		nsec = nsec%86400;
		double fsec = (stt_smjd%2)*0.5+stt_offs;
		nsec += (long int)fsec;
		fsec = fsec-(long int)fsec;

		MJD tt;
		tt.stt_imjd = nimjd;
		tt.stt_smjd = nsec;
		tt.stt_offs = fsec;
		return tt;
	}
	const MJD & operator =(const MJD &t)
	{
		stt_imjd = t.stt_imjd;
		stt_smjd = t.stt_smjd;
		stt_offs = t.stt_offs;
		return *this;
	}
	const MJD & operator =(long double mjd)
	{
		format(mjd);
		return *this;
	}

	bool operator <(const MJD &t)
	{
		if (stt_imjd < t.stt_imjd)
		{
			return true;
		}
		else if (stt_imjd == t.stt_imjd)
		{
			if (stt_smjd < t.stt_smjd)
			{
				return true;
			}
			else if (stt_smjd == t.stt_smjd)
			{
				if (stt_offs < t.stt_offs)
				{
					return true;
				}
			}
		}

		return false;
	}

	bool operator >(const MJD &t)
	{
		if (stt_imjd > t.stt_imjd)
		{
			return true;
		}
		else if (stt_imjd == t.stt_imjd)
		{
			if (stt_smjd > t.stt_smjd)
			{
				return true;
			}
			else if (stt_smjd == t.stt_smjd)
			{
				if (stt_offs > t.stt_offs)
				{
					return true;
				}
			}
		}

		return false;
	}

	bool operator ==(const MJD &t)
	{
		if (stt_imjd == t.stt_imjd and stt_smjd == t.stt_smjd and stt_offs == t.stt_offs)
			return true;
		else
			return false;
	}

	MJD operator +(const MJD &t)
	{
		MJD tt;
		tt.stt_imjd = stt_imjd + t.stt_imjd;
		tt.stt_smjd = stt_smjd + t.stt_smjd;
		tt.stt_offs = stt_offs + t.stt_offs;
		return tt;
	}
	MJD operator +(double sec)
	{
		long int nsec = floor(sec);
		double fsec = sec-nsec;
		MJD tt;
		tt.stt_imjd = stt_imjd;
		tt.stt_smjd = stt_smjd + nsec;
		tt.stt_offs = stt_offs + fsec;
		return tt;
	}
	void operator +=(double sec)
	{
		long int nsec = floor(sec);
		double fsec = sec-nsec;
		stt_smjd += nsec;
		stt_offs += fsec;
	}
	MJD operator -(const MJD &t)
	{
		MJD tt;
		tt.stt_imjd = stt_imjd - t.stt_imjd;
		tt.stt_smjd = stt_smjd - t.stt_smjd;
		tt.stt_offs = stt_offs - t.stt_offs;
		return tt;
	}
	MJD operator -(double sec)
	{
		long int nsec = floor(sec);
		double fsec = sec-nsec;
		MJD tt;
		tt.stt_imjd = stt_imjd;
		tt.stt_smjd = stt_smjd - nsec;
		tt.stt_offs = stt_offs - fsec;
		return tt;
	}
	void operator -=(double sec)
	{
		long int nsec = floor(sec);
		double fsec = sec-nsec;
		stt_smjd -= nsec;
		stt_offs -= fsec;
	}

	bool gregorian(struct tm* gregdate, double* fsec) const
	{
		int julian_day = stt_imjd + 2400001;

		int n_four = 4  * (julian_day+((6*((4*julian_day-17918)/146097))/4+1)/2-37);
		int n_dten = 10 * (((n_four-237)%1461)/4) + 5;

		gregdate->tm_year = n_four/1461 - 4712 - 1900; // extra -1900 for C struct tm
		gregdate->tm_mon  = (n_dten/306+2)%12;         // struct tm mon 0->11
		gregdate->tm_mday = (n_dten%306)/10 + 1;

		ss2hhmmss(&gregdate->tm_hour, &gregdate->tm_min, &gregdate->tm_sec, stt_smjd);

		if (fsec)
		*fsec = stt_offs;

		gregdate->tm_isdst = -1;
		time_t date = mktime (gregdate);
		if (date == (time_t)-1)
		return false;

		return true;
	}
public:
	static void ss2hhmmss(int* hours, int* min, int* sec, int seconds)
	{
	  *hours   = seconds/3600;
	  seconds -= 3600 * (*hours);
	  *min     = seconds/60;
	  seconds -= 60 * (*min);
	  *sec     = seconds;
	}

	friend bool operator <(const MJD &t1, const MJD &t2)
	{
		if (t1.stt_imjd < t2.stt_imjd)
		{
			return true;
		}
		else if (t1.stt_imjd == t2.stt_imjd)
		{
			if (t1.stt_smjd < t2.stt_smjd)
			{
				return true;
			}
			else if (t1.stt_smjd == t2.stt_smjd)
			{
				if (t1.stt_offs < t2.stt_offs)
				{
					return true;
				}
			}
		}

		return false;
	}

	friend bool operator >(const MJD &t1, const MJD &t2)
	{
		if (t1.stt_imjd > t2.stt_imjd)
		{
			return true;
		}
		else if (t1.stt_imjd == t2.stt_imjd)
		{
			if (t1.stt_smjd > t2.stt_smjd)
			{
				return true;
			}
			else if (t1.stt_smjd == t2.stt_smjd)
			{
				if (t1.stt_offs > t2.stt_offs)
				{
					return true;
				}
			}
		}

		return false;
	}

public:
	long int stt_imjd;
	long int stt_smjd;
	double stt_offs;
};

#endif /* MJD_H_ */
