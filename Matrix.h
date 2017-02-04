/*		A 3x3 matrix package for C++
		(c) Walter Brisken 5-20-1993

			Updated 6-28-1996       Fremont, CA
                        Updated 11-4-1996       Princeton, NJ
                        Updated 8-16-1997       Princeton, NJ
                        Updated 12-10-1998      Socorro, NM     (templatized)

*/

#ifndef __MATRIX___
#define __MATRIX___

#include "Vector.h"

template <class T> class Matrix
{
public:
	T a[3][3];  // (down, across)

	Matrix(const T &a11 = 0, const T &a12 = 0, const T &a13 = 0,
	       const T &a21 = 0, const T &a22 = 0, const T &a23 = 0,
	       const T &a31 = 0, const T &a32 = 0, const T &a33 = 0) 
	{
		a[0][0] = a11; 
		a[1][0] = a21; 
		a[2][0] = a31;
		a[0][1] = a12; 
		a[1][1] = a22; 
		a[2][1] = a32;
		a[0][2] = a13; 
		a[1][2] = a23; 
		a[2][2] = a33;
	}

	Matrix(const Vector<T> &v1, const Vector<T> &v2, const Vector<T> &v3)
	{
		a[0][0] = v1.x;
		a[1][0] = v1.y;
		a[2][0] = v1.z;
		a[0][1] = v2.x;
		a[1][1] = v2.y;
		a[2][1] = v2.z;
		a[0][2] = v3.x;
		a[1][2] = v3.y;
		a[2][2] = v3.z;
	}

	Matrix(const Vector<T> &v)
	{
		a[0][0] = v.x;
		a[1][0] = 0;
		a[2][0] = 0;
		a[0][1] = 0;
		a[1][1] = v.y; 
		a[2][1] = 0; 
		a[0][2] = 0;  
		a[1][2] = 0;
		a[2][2] = v.z;
	}

	Matrix(const Matrix &arg)
	{
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++) 
				a[i][j] = arg.a[i][j];
	}

	Matrix(double rotx, double roty, double rotz)
	{
		double sx, cx, sy, cy, sz, cz;
		sx = sin(rotx);
		cx = cos(rotx);
		sy = sin(roty);
		cy = cos(roty);
		sz = sin(rotz);
		cz = cos(rotz);
		a[0][0] = cz*cy;
		a[0][1] = -sz*cx+cz*sy*sx;
		a[0][2] = sz*sx+cz*sy*cx;
		a[1][0] = sz*cy;
		a[1][1] = cz*cx+sz*sy*sx;
		a[1][2] = -cz*sx+sz*sy*cx;
		a[2][0] = -sy;
		a[2][1] = cy*sx;
		a[2][2] = cy*cx;
	}


	friend Matrix operator + (const Matrix &m, const Matrix &arg)
	{
		Matrix result;
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++) 
				result.a[i][j] = m.a[i][j] + arg.a[i][j];
		return result;
	}
	
	friend Matrix operator - (const Matrix &m, const Matrix &arg)
	{
		Matrix result;
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++) 
				result.a[i][j] = m.a[i][j] - arg.a[i][j];
		return result;
	}
	
	friend Matrix operator - (const Matrix &m)
	{
		Matrix result;
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++) 
				result.a[i][j] = -m.a[i][j];
		return result;
	}
	
	Matrix &operator = (const Matrix &arg)
	{
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++)
				a[i][j] = arg.a[i][j];
		return *this;
	}

	friend Matrix operator * (const Matrix &m, const Matrix &arg)
	{
		Matrix res;
		for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				res.a[i][j] += m.a[i][k]*arg.a[k][j];
		return res;
	}

	friend Matrix operator * (const Matrix &m, T arg)
	{
		Matrix res;
		for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++)
			res.a[i][j] = arg*m.a[i][j];
		return res;
	}

	friend Vector<T> operator * (const Matrix &m, const Vector<T> &arg)
	{
		Vector<T> res;
		res.x = m.a[0][0]*arg.x + m.a[0][1]*arg.y + m.a[0][2]*arg.z;
		res.y = m.a[1][0]*arg.x + m.a[1][1]*arg.y + m.a[1][2]*arg.z;
		res.z = m.a[2][0]*arg.x + m.a[2][1]*arg.y + m.a[2][2]*arg.z;
		return res;
	}

	friend Vector<T> operator * (const Vector<T> &arg, const Matrix &m)
	{
		Vector<T> res;
		res.x = m.a[0][0]*arg.x + m.a[1][0]*arg.y + m.a[2][0]*arg.z;
		res.y = m.a[0][1]*arg.x + m.a[1][1]*arg.y + m.a[2][1]*arg.z;
		res.z = m.a[0][2]*arg.x + m.a[1][2]*arg.y + m.a[2][2]*arg.z;
		return res;
	}

//	friend Matrix operator ^ (const Matrix &m, int power);

	friend Matrix operator ~ (const Matrix &m)     // Transpose
	{
		Matrix res;
		for(int i = 0; i < 3; i++) 
			for(int j = 0; j < 3; j++) 
				res.a[i][j] = m.a[j][i];
		return res;
	}

	T det()
	{
		return  a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+
                	a[0][2]*a[1][0]*a[2][1]-a[2][0]*a[1][1]*a[0][2]-
                	a[2][1]*a[1][2]*a[0][0]-a[2][2]*a[1][0]*a[0][1];
	}

	T trace()
	{
		return a[0][0] + a[1][1] + a[2][2];
	}

	void symmetrize()
	{
		T sum;
		sum = a[0][1] + a[1][0];
		a[0][1] = a[1][0] = sum/2.0;
		sum = a[1][2] + a[2][1];
		a[1][2] = a[2][1] = sum/2.0;
		sum = a[2][0] + a[0][2];
		a[2][0] = a[0][2] = sum/2.0;
	}

	Matrix adjoint()
	{
		return   Matrix(a[1][1]*a[2][2]-a[1][2]*a[2][1],
                        	a[0][2]*a[2][1]-a[0][1]*a[2][2],
                        	a[0][1]*a[1][2]-a[0][2]*a[1][1],
                        	a[1][2]*a[2][0]-a[1][0]*a[2][2],
                        	a[0][0]*a[2][2]-a[0][2]*a[2][0],
                        	a[0][2]*a[1][0]-a[0][0]*a[1][2],
                        	a[1][0]*a[2][1]-a[1][1]*a[2][0],
                        	a[0][1]*a[2][0]-a[0][0]*a[2][1],
                        	a[0][0]*a[1][1]-a[0][1]*a[1][0]);

	}

	Matrix inv()
	{
		T d = det();
		if(d == 0.0) return Matrix();
		return adjoint()*(1.0/d);
	}

	friend Matrix outerproduct(const Vector<T> &v1, const Vector<T> &v2)
	{
		return Matrix(v1.x*v2.x, v1.x*v2.y, v1.x*v2.z,
			      v1.y*v2.x, v1.y*v2.y, v1.y*v2.z,
			      v1.z*v2.x, v1.z*v2.y, v1.z*v2.z);
	}

	friend ostream &operator << (ostream &s, const Matrix &m)
	{
		return s << "[ [" << m.a[0][0] << ", "   << m.a[0][1] << ", " 
				  << m.a[0][2] << "], [" << m.a[1][0] << ", "
				  << m.a[1][1] << ", "   << m.a[1][2] << "], ["
				  << m.a[2][0] << ", "   << m.a[2][1] << ", "
				  << m.a[2][2] << "] ]";
	}
};

#endif
