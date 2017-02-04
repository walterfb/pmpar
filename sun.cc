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

#include "sun.h"
#include "coordutils.h"

/* Algorithm from 2003 Astronomical almanac.
 * answer is returned as a 3 vector for convenience
 */

double quickmod(double a, double b)
{
	while(a < 0.0)
	{
		a += b;
	}
	while(a >= b)
	{
		a -= b;
	}

	return a;
}

// Returns vector in direction of sun for given MJD.

Vector<double> sunpos(double day)
{
	double n, L, grad;
	double lambdarad, epsrad;
	double alpharad, deltarad;
	double audist;
	
	n = day - 51544.5;
	L = quickmod(280.460 + 0.9856474*n, 360.0);
	grad = quickmod(357.528 + 0.9856003*n, 360.0)*PI_180;

	lambdarad = quickmod(L + 1.915*sin(grad) + 0.020*sin(2.0*grad), 360.0)*PI_180;
	epsrad = (23.439 - 0.0000004*n)*PI_180;
	alpharad = atan(cos(epsrad)*tan(lambdarad));
	while(alpharad - lambdarad > PI_2) alpharad -= M_PI;	// tan ambig
	while(alpharad - lambdarad < -PI_2) alpharad += M_PI;
	deltarad = asin(sin(epsrad)*sin(lambdarad));

	audist = 1.00014 - 0.01671*cos(grad) - 0.00014*cos(2.0*grad);

	return radec2vec(alpharad, deltarad);
}
