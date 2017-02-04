/*  Copyright (C) 1993-2017 Walter Brisken  wbrisken@lbo.us
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

/*			A Vector Package for C++
	  	    (C)  Walter Brisken 1-17-1993
			    Pasadena, CA

			Updated 6-28-96
			Updated 11-4-96
			Updated 8-30-97
			Updated 12-7-98	Templatized
*/

#include <cmath>
#include "Vector.h"

template<class T>
Vector<T> Vector<T>::rotate(double cos1, double sin1, double cos2, double sin2)
{
	Vector<T> result;
	double temp;
	temp = x*sin1 + z*cos1;
	result.x =  x*cos1 - z*sin1;
	result.y = -y*cos2 + temp*sin2;
	result.z = -y*sin2 - temp*cos2;
	return result;
}

template<class T>
Vector<T> &Vector<T>::rotatex(double cos1, double sin1)
{
	double temp = y;
	y = temp*cos1 - z*sin1;
	z = temp*sin1 + z*cos1;
	return *this;
}

template<class T>
Vector<T> &Vector<T>::rotatey(double cos1, double sin1)
{
	double temp = x;
	x = temp*cos1 - z*sin1;
	z = temp*sin1 + z*cos1;
	return *this;
}

template<class T>
Vector<T> &Vector<T>::rotatez(double cos1, double sin1)
{
	double temp = x;
	x = temp*cos1 - y*sin1;
	y = temp*sin1 + y*cos1;
	return *this;
}

template<class T>
Vector<T> Vector<T>::rev_rotate(double cos1, double sin1, 
				double cos2, double sin2)
{
	Vector<T> result;
	double temp;
	temp = y*sin2 + z*cos2;
	result.y = y*cos2 - z*sin2;
	result.x = -x*cos1 + temp*sin1;
	result.z = -x*sin1 - temp*cos1;
	return result;
}
