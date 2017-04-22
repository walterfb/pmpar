/*  Copyright (C) 2000-2017 Walter Brisken  wbrisken@lbo.us
 *
 *  This file is part of pmpar.
 *
 *  pmpar is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  pmpar is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with pmpar.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cstdio>
#include "coordutils.h"

Vector<double> radec2vec(double alpharad, double deltarad)
{
	double cd = cos(deltarad);

	return Vector<double>(cos(alpharad)*cd, sin(alpharad)*cd, sin(deltarad));
}

// Assumes J2000 equinox, all coords in radians

void equ2gal(double ra, double dec, double *lout, double *bout)
{
	const double stheta = 0.88998808748;
	const double ctheta = 0.45598377618;
	const double psi = 0.57477043300;
	const double phi = 4.9368292465;

	double a, b, sb, cb, cbsa;

	a = ra - phi;
	b = dec;
	sb = sin(b);
	cb = cos(b);
	cbsa = cb * sin(a);
	b = -stheta * cbsa + ctheta * sb;
	*lout = atan2( ctheta * cbsa + stheta * sb, cb * cos(a))+psi;
	*bout = asin(b);
	while(*lout < 0.0)
	{
		*lout = *lout + 2.0*M_PI;
	}
	while(*lout >= 2.0*M_PI)
	{
		*lout = *lout - 2.0*M_PI;
	}
}

void muequ2gal(double ra, double dec, double mu_a, double mu_d, double *mu_l, double *mu_b)
{
	double l0, b0, l, b;
	double eps = 0.0001;
	double l_a, b_a, l_d, b_d;

	equ2gal(ra, dec, &l0, &b0);
	equ2gal(ra+eps/cos(dec), dec, &l, &b);
	l_a = (l - l0)*cos(b0) / eps;
	b_a = (b - b0) / eps;
	
	// The Richard Dodson fix: l0, b0 -> l, b below
	equ2gal(ra, dec+eps, &l, &b);
	
	l_d = (l - l0)*cos(b0) / eps;
	b_d = (b - b0) / eps;
	*mu_l = mu_a*l_a + mu_d*l_d;
	*mu_b = mu_a*b_a + mu_d*b_d;
}

double text2double(const char *str)
{
	int colons = 0;
	int len;
	double sign = 1.0;
	char str2[100], *str3;
	double h = 0, m, s;
	double ret = -999999.999;
	
	for(len = 0; (str2[len] = str[len]) != 0; ++len)
	{
		if(str[len] == ':') ++colons;
	}

	if(str[0] == '-' || str[0] == '+')
	{
		if(str[0] == '-')
		{
			sign = -1.0;
		}
		str3 = str2+1;
	}
	else
	{
		str3 = str2;
	}

	if(colons == 0) 
	{
#if 0
		int digs;
		for(digs = 0; str3[digs] != 0 && str3[digs] != '.'; digs++);
		switch(digs)
		{
			case 0:
			case 1:
			case 2:
				if(sscanf(str, "%lf", &h) == 1) ret = h;
				break;
			case 4:
				h = 10.0*(str3[0] - '0') + (str3[1] - '0');
				if(sscanf(str3 + 2, "%lf", &m) == 1)
					ret = sign*(h + m/60.0);
				break;
			case 7:
				h = 100.0*(str3[0] - '0');
				str3++;
			case 6:
				h += 10.0*(str3[0] - '0') + (str3[1] - '0');
				m =  10.0*(str3[2] - '0') + (str3[3] - '0');
				if(sscanf(str3+4, "%lf", &s) == 1)
					ret = sign*(h + m/60.0 + s/3600.0);
				break;
		}
#endif
		sscanf(str, "%lf", &h);
		ret = h/3600.0;
	}
	else if(colons == 2)
	{
		if(sscanf(str, "%lf:%lf:%lf", &h, &m, &s) == 3) 
		{
			ret = fabs(h) + m/60.0 + s/3600.0;
			if(str[0] == '-')
			{
				ret = -ret;
			}
		}
	}

	return ret;
}

void double2text(double num, char *str, char sep)
{
	int h, m;

	if(num < 0)
	{
		num = -num;
		str[0] = '-';
		str++;
	}
	h = (int)num;
	num -= h;
	num *= 60.0;
	m = (int)num;
	num -= m;
	num *= 60.0;
	sprintf(str, "%02d%c%02d%c%09.6f", h, sep, m, sep, num);
}
