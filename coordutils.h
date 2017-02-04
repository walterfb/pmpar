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

#ifndef __COORDUTILS_H__
#define __COORDUTILS_H__

#include "Vector.h"
#include "math.h"

#define PI_180  (M_PI/180.0)
#define PI_12   (M_PI/12.0)
#define PI2     (M_PI*2.0)
#define PI_2    (M_PI/2.0)
#define MAS_RAD	(3600000.0*180.0/M_PI) 

void equ2gal(double ra, double dec, double *lout, double *bout);

void muequ2gal(double ra, double dec, double mu_a, double mu_d, double *mu_l, double *mu_b);	

Vector<double> radec2vec(double alpharad, double deltarad);

double text2double(const char *str);

void double2text(double num, char *str, char sep = ':');

#endif
