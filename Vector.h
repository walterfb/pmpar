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
		      (c)  Walter Brisken 1-17-1993
				Pasadena, CA
				Princeton, NJ
				Socorro, NM
				St Paul, MN
*/

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <iostream>
#include <cmath>

template <class T> class Vector
{
public:
	T x, y, z;

	Vector(T x1 = 0, T y1 = 0, T z1 = 0) : x(x1), y(y1), z(z1) { }
	Vector(const Vector &v) : x(v.x), y(v.y), z(v.z) { }

	Vector &operator = (const Vector &v)	// equation
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	Vector &operator += (const Vector &v)	// a=a+b
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	Vector &operator -= (const Vector &v)	// a=a-b
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	Vector &operator *= (const T &a)		// a=a*b
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	Vector &operator *= (const Vector &v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}

	Vector &operator /= (const T &a)		// a=a/b
	{
		x /= a;
		y /= a;
		z /= a;
		return *this;
	}

	friend Vector operator - (const Vector &v)      	// negation
	{
		return Vector<T>(-v.x, -v.y, -v.z);
	}
	
	friend Vector operator + (const Vector &a, const Vector &b)
	{
		return Vector<T>(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	
	friend Vector operator - (const Vector &a, const Vector &b)
	{
		return Vector<T>(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	
	friend Vector operator * (const Vector &a, const Vector &b)
	{
		return Vector<T>(a.x*b.x, a.y*b.y, a.z*b.z);
	}
	
	friend Vector operator * (const Vector &v, const T &a) 
	{
		return Vector<T>(v.x*a, v.y*a, v.z*a);
	}

	friend Vector operator / (const Vector &a, const Vector &b)
	{
		return Vector(a.x/b.x, a.y/b.y, a.z/b.z);
	}

	friend Vector operator / (const Vector &v, T a)
	{
		return Vector(v.x/a, v.y/a, v.z/a);
	}

	// dot product
	friend T operator % (const Vector &u, const Vector &v)
	{
		return u.x*v.x + u.y*v.y + u.z*v.z;
	}

	// cross product
	friend Vector operator ^ (const Vector &u, const Vector &v)
	{
		return Vector<T>(u.y*v.z - u.z*v.y, 
				 u.z*v.x - u.x*v.z, 
				 u.x*v.y - u.y*v.x);
	}

	friend T Dot(const Vector &u, const Vector &v)
	{
		return u.x*v.x + u.y*v.y + u.z*v.z;
	}

	friend T Square(const Vector &v)
	{
		return v.x*v.x + v.y*v.y + v.z*v.z;
	}

	T length()
	{
		return sqrt(x*x + y*y + z*z);
	}

	Vector &normalize()
	{
		T l = length();
		if(l < 0.000001) l=0.000001;
		x/=l;
		y/=l;
		z/=l;
		return *this;
	}

	Vector proj(const Vector &v)
	{
		T div = Square(v);
		T fac;
		if(div < 0.0000001) div = 0.0000001;
		fac = (x*v.x + y*v.y + z*v.z)/div;
		return Vector<T>(v.x*fac,  v.y*fac, v.z*fac);
	}

	Vector unitproj(const Vector &unitv)
	{
		T fac = x*unitv.x + y*unitv.y + z*unitv.z;
		return Vector<T>(unitv.x*fac, unitv.y*fac, unitv.z*fac);
	}

	Vector ortho(const Vector &v)
	{
		return *this-proj(v);
	}

	inline Vector operator~()            	// normalation
	{
		T l = sqrt(x*x + y*y + z*z);
		if(l < 0.000001) return Vector(1, 0, 0);
		else return Vector<T>(x/l, y/l, z/l);
	}

//	friend ostream& operator << (ostream &s, const Vector &v)
//	{
//		return s << "(" << v.x << ", " << v.y << ", " << v.z << ")";
//	}

	Vector min(const Vector &);		// and much, much more!!
	Vector max(const Vector &);
	Vector rotate(double cos1, double sin1, double cos2, double sin2);
	Vector &rotatex(double cos1, double sin1);
	Vector &rotatey(double cos1, double sin1);
	Vector &rotatez(double cos1, double sin1);
	Vector rev_rotate(double cos1, double sin1, double cos2, double sin2);
};

#endif
