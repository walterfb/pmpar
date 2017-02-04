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
#include <cstdlib>
#include <cstring>
#include "sun.h"
#include "coordutils.h"
#include "config.h"

#define MAXEPOCHS	400

void getab(double day, double ra, double dec, double *pi_a, double *pi_d);
void printday(double mjd, int opt = 0);

int epochs = 0;
double epoch[MAXEPOCHS];
double a[MAXEPOCHS];
double d[MAXEPOCHS];
double siga[MAXEPOCHS];
double sigd[MAXEPOCHS];
int inplace = 0;

char parmnames[9][9] =
{
	"  RA    ", 
	"  Dec   ", 
	"  mu_a  ", 
	"  mu_d  ", 
	"  pi    ", 
	"  sh_x  ", 
	"  sh_y  ", 
	"  dmu_x ", 
	"  dmu_y "
};

/* Move all internals to use true milliarcsecond coordinates.  Used to have
mas/cos(dec) for RA. */

class Pulsar
{
public:
	double a_0, d_0, siga_0, sigd_0;
	double mu_a, mu_d, sigmu_a, sigmu_d;
	double sh_x, sh_y, sig_sh_x, sig_sh_y;
	double dmu_x, dmu_y, sig_dmu_x, sig_dmu_y;
	double pi, sigpi, dm;
	double reverse;
	double epoch;
	double breaktime;
	double normcov[9][9];
	int fit[9];
	double parm[9];
	double sigmu_a_given;
	double sigmu_d_given;
	double sigpi_given;
	char name[100], ref[100];

	Pulsar() 
	{ 
		dm = -1.0;
		a_0 = d_0 = mu_a = mu_d = pi = 0.0; 
		siga_0 = sigd_0 = sigmu_a = sigmu_d = sigpi = 0.0;
		sh_x = sh_y = sig_sh_x = sig_sh_y = 0.0;
		epoch = 2000.0;
		breaktime = -1.0;
		for(int i = 0; i < 9; ++i)
		{
			fit[i] = 1;
			parm[i] = 0.0;
		}
		name[0] = 0;
		ref[0] = 0;
		sigmu_a_given = 0.0;
		sigmu_d_given = 0.0;
		sigpi_given = 0.0;
		reverse = 1.0;
	}	

	Pulsar &operator = (const Pulsar &p)
	{
		a_0     = p.a_0;
		d_0     = p.d_0;
		mu_a    = p.mu_a;
		mu_d    = p.mu_d;
		pi      = p.pi;
		sh_x 	= p.sh_x;
		sh_y 	= p.sh_y;
		dmu_x	= p.dmu_x;
		dmu_y	= p.dmu_y;
		epoch   = p.epoch;
		siga_0  = p.siga_0;
		sigd_0  = p.sigd_0;
		sigmu_a = p.sigmu_a;
		sigmu_d = p.sigmu_d;
		sigpi   = p.sigpi;
		sig_sh_x= p.sig_sh_y;
		sig_sh_y= p.sig_sh_y;
		sig_dmu_x=p.sig_dmu_x;
		sig_dmu_y=p.sig_dmu_y;
		sigmu_a_given = p.sigmu_a_given;
		sigmu_d_given = p.sigmu_d_given;
		sigpi_given = p.sigpi_given;
		reverse = p.reverse;

		for(int i = 0; i < 9; ++i)
		{
			fit[i] = p.fit[i];
		}
		for(int i = 0; i < 9; ++i)
		{
			parm[i] = p.parm[i];
		}
		for(int i = 0; i < 9; ++i)
		{
			for(int j = 0; j < 9; ++j)
			{
				normcov[i][j] = p.normcov[i][j];
			}
		}

		return *this;
	}

	// Barycenter coords as func of t, in radians
	double ra_sun(double t)  
	{ 
		double res;
		
		res = a_0 + mu_a*(t-epoch)/cos(d_0)/365.242175; 
		if(epoch < breaktime)
		{
			res += (sh_x + dmu_x*(t-epoch))/cos(d_0);
		}

		return res;
	}
	double dec_sun(double t) 
	{ 
		double res;
		
		res = d_0 + mu_d*(t-epoch)/365.242175; 
		if(epoch < breaktime)
		{
			res += (sh_y + dmu_y*(t-epoch));
		}

		return res;
	}
	double ra_earth(double t) 
	{
		double ras, decs, pi_a, pi_d;
		
		ras = ra_sun(t);
		decs = dec_sun(t);
		getab(t, ras, decs, &pi_a, &pi_d);
		
		return ras + pi*pi_a/cos(d_0);
	}
	double dec_earth(double t)
	{
		double ras, decs, pi_a, pi_d;
		
		ras = ra_sun(t);
		decs = dec_sun(t);
		getab(t, ras, decs, &pi_a, &pi_d);
		
		return decs + pi*pi_d;
	}
	void maspos(double t, int dropm, int dropp, int drops, double *x, double *y)
	{
		double a1 = 0.0, d1 = 0.0;
		
		if(t < 0.0)
		{
			t = epoch;
		}
		if(!dropm)
		{
			a1 += mu_a*(t-epoch)/365.242175;
			d1 += mu_d*(t-epoch)/365.242175;
		}
		if(!dropp)
		{
			double pi_a, pi_d;
			double ras = ra_sun(t);
			double decs = dec_sun(t);

			getab(t, ras, decs, &pi_a, &pi_d);
			a1 += pi*pi_a;
			d1 += pi*pi_d;
		}
		if(!drops) if(t < breaktime)
		{
			a1 += (sh_x + dmu_x*(t-epoch)/365.242175);
			d1 += (sh_y + dmu_y*(t-epoch)/365.242175);
		}
		*x = a1*MAS_RAD;
		*y = d1*MAS_RAD;
	}
	void maspoint(int e, int dropm, int dropp, int drops, double *x, double *y)
	{
		double a1 = (a[e]-a_0)*cos(d_0), d1 = d[e]-d_0;
		double t = ::epoch[e];

		if(dropm)
		{
			a1 -= mu_a*(t-epoch)/365.242175;
			d1 -= mu_d*(t-epoch)/365.242175;
		}
		if(dropp)
		{
			double pi_a, pi_d;
			double ras = ra_sun(t);
			double decs = dec_sun(t);

			getab(t, ras, decs, &pi_a, &pi_d);
			a1 -= pi*pi_a;
			d1 -= pi*pi_d;
		}
		if(drops) if(t < breaktime)
		{
			a1 += (sh_x + dmu_x*(t-epoch)/365.242175);
			d1 += (sh_y + dmu_y*(t-epoch)/365.242175);
		}
		*x = a1*MAS_RAD;
		*y = d1*MAS_RAD;
	}

	void printnormcov()
	{
		printf("Normalized covariance matrix:\n");
		printf("        ");
		for(int i = 0; i < 9; ++i)
		{
			if(fit[i])
			{
				printf("%s", parmnames[i]);
			}
		}
		printf("\n");
		for(int j = 0; j < 9; ++j)
		{
			if(fit[j])
			{
				printf("%s", parmnames[j]);
				for(int i = 0; i < 9; ++i)
				{
					if(fit[i])
					{
						printf("%7.4f ", normcov[i][j]);
					}
				}
				printf("\n");
			}
		}
		printf("\n");
	}
	
	void printstatistics()
	{
		double dx, dy, x, y, mx, my, X2=0.0, sx=0.0, sy=0.0;
		int dof = 0;

		for(int i = 0; i < 9; ++i)
		{
			if(fit[i] == 1)
			{
				dof++;
			}
		}

		for(int e = 0; e < epochs; e++)
		{
			maspos(e, 1, 1, 1, &mx, &my);
			maspoint(e, 1, 1, 1, &x, &y);
			dx = x - mx;
			dy = y - my;
			sx += dx*dx;
			sy += dy*dy;
			X2 += (dx*dx)/(siga[e]*siga[e]*cos(d[e])*cos(d[e])*MAS_RAD*MAS_RAD);
			X2 += (dy*dy)/(sigd[e]*sigd[e]*MAS_RAD*MAS_RAD);
		}

		if(dof > 0)
		{
			printf("\nNDOF = %d\n", dof);
			printf("Scatterx = %f mas\n", sqrt(sx/(double)epochs));
			printf("Scattery = %f mas\n", sqrt(sy/(double)epochs));
			printf("Reduced Chi^2 = %f\n", X2/(double)(2*epochs-dof));
		}
	}

	void saveaspsr()
	{
		double a, b;
		char str[100];
		double mu_l, mu_b;
		FILE *out;
	
		sprintf(str, "%s.psr", name);
		out = fopen(str, "w");
		if(!out)
		{
			fprintf(stderr, "Cannot open %s for write\n", str);
			return;
		}
		fprintf(out, "# This file generated by pmpar\n");
		
		fprintf(out, "Name = %s\n", name);
		fprintf(out, "Epoch = %9.3f\n", epoch);
		fprintf(out, "ra = %f\n", a_0);
		fprintf(out, "dec = %f\n", d_0);
		equ2gal(a_0, d_0, &a, &b);
		muequ2gal(a_0, d_0, mu_a, mu_d, &mu_l, &mu_b);
		fprintf(out, "gall = %f\n", a);
		fprintf(out, "galb = %f\n", b);
		fprintf(out, "mu_a = %f # mu_ra*cos(dec)\n", mu_a*MAS_RAD);
		fprintf(out, "mu_d = %f\n", mu_d*MAS_RAD);
		if(sigmu_a > 0.0) 
		{
			fprintf(out, "sigmu_a = %f\n", sigmu_a*MAS_RAD);
		}
		else
		{
			fprintf(out, "sigmu_a = %f\n", sigmu_a_given);
		}
		if(sigmu_d > 0.0)
		{
			fprintf(out, "sigmu_d = %f\n", sigmu_d*MAS_RAD);
		}
		else
		{
			fprintf(out, "sigmu_d = %f\n", sigmu_d_given);
		}
		fprintf(out, "pi = %f\n", pi*MAS_RAD);
		if(sigpi > 0.0) 
		{
			fprintf(out, "sigpi = %f\n", sigpi*MAS_RAD);
		}
		else
		{
			fprintf(out, "sigpi = %f\n", sigpi_given);
		}
		fprintf(out, "DM = %f\n", dm);
		fclose(out);
	}

	void print()
	{
		double a, b;
		char str[100];
		double dist0, distp, distm, mu_l, mu_b;
	
		if(name[0] != 0)
		{
			printf("Name  = %s\n", name);
		}
		if(ref[0] != 0)
		{
			printf("Ref.  = %s\n", ref);
		}
		printf("epoch = %9.3f = ", epoch);
		printday(epoch);
		double2text(a_0/PI_12, str);
		printf("RA    = %s +- %f\n", str, siga_0*3600.0/PI_12);
		double2text(d_0/PI_180, str);
		printf("Dec   = %s +- %f\n", str, sigd_0*3600.0/PI_180);
		equ2gal(a_0, d_0, &a, &b);
		muequ2gal(a_0, d_0, mu_a, mu_d, &mu_l, &mu_b);
		printf("l     = %f degrees\n", a/PI_180);
		printf("b     = %f degrees\n", b/PI_180);
		printf("mu_a  = %f +- %f mas/yr  (mu_ra*cos(dec))\n", mu_a*MAS_RAD, sigmu_a*MAS_RAD);
		printf("mu_d  = %f +- %f mas/yr\n", mu_d*MAS_RAD, sigmu_d*MAS_RAD);
		printf("mu_l  = %f mas/yr\n", mu_l*MAS_RAD);
		printf("mu_b  = %f mas/yr\n", mu_b*MAS_RAD);
		printf("pi    = %f +- %f mas\n", pi*MAS_RAD, sigpi*MAS_RAD);
		dist0 = 1000.0/(pi*MAS_RAD);
		distp = 1000.0/((pi-sigpi)*MAS_RAD);
		distm = 1000.0/((pi+sigpi)*MAS_RAD);
		printf("dist  = %f + %f - %f pc\n", dist0, distp-dist0, dist0-distm);
		if(dm >= 0.0)
		{
			printf("ne    = %f +- %f cm^-3\n", dm/dist0, (dm/distm - dm/distp)/2.0);
		}
		printf("v_t   = %f km/sec\n", 0.00474*dist0*MAS_RAD*sqrt(mu_a*mu_a+mu_d*mu_d));
		if(breaktime > 0.0)
		{
			printf("shift x = %f +- %f mas\n", sh_x*MAS_RAD, sig_sh_x*MAS_RAD);
			printf("shift y = %f +- %f mas\n", sh_y*MAS_RAD, sig_sh_y*MAS_RAD);
			printf("dmu_a = %f +- %f mas/yr\n", dmu_x*MAS_RAD, sig_dmu_x*MAS_RAD);
			printf("dmu_d = %f +- %f mas/yr\n", dmu_y*MAS_RAD, sig_dmu_y*MAS_RAD);
		}
		printstatistics();
		printf("\n");
	}

	void printminimal()
	{
		printf("%10.8f %10.8f  %10.8f %10.8f  %f %f  %f %f  %f %f\n",
			a_0*3600.0/PI_12, siga_0*3600.0/PI_12,	  /* RA and error in seconds of time */
			d_0*3600.0/PI_180, sigd_0*3600.0/PI_180,  /* Dec and error in arcsec */
			mu_a*MAS_RAD, sigmu_a*MAS_RAD,		  /* mu_a in real arcsec and error */
			mu_d*MAS_RAD, sigmu_d*MAS_RAD,		  /* mu_d and error in arcsec */
			pi*MAS_RAD, sigpi*MAS_RAD);
	}
};

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

void gaussj9(double a[9][9], double b[9])
{
	int indxc[9], indxr[9], ipiv[9];
        int icol=0, irow=0;
        double big, dum, pivinv;

        for(int j = 0; j < 9; ++j)
	{
		ipiv[j] = 0;
	}
        for(int i = 0; i < 9; ++i)
        {
                big = 0.0;
                for(int j = 0; j < 9; ++j) 
		{
			if(ipiv[j] != 1) 
			{
				for(int k = 0; k < 9; ++k)
				{
				        if(ipiv[k] == 0)
				        {
				                if(fabs(a[j][k]) >= big)
				                {
				                        big = fabs(a[j][k]);
						        irow=j;
							icol=k;
						}
					}
					else if(ipiv[k] > 1)
					{
						return;
					}
				}
				ipiv[icol]++;
				if(irow != icol)
				{
					for(int l = 0; l < 9; ++l)
					{
						SWAP(a[irow][l],a[icol][l])
					}
					SWAP(b[irow], b[icol])
				}
			        indxr[i]=irow;
				indxc[i]=icol;
			        if(a[icol][icol] == 0.0)
				{
					return;
				}
				pivinv = 1.0/a[icol][icol];
			        a[icol][icol] = 1.0;
			        for(int l = 0; l < 9; l++)
				{
					a[icol][l] *= pivinv;
				}
			        b[icol] *= pivinv;
				for(int ll = 0; ll < 9; ++ll)
				{
					if(ll != icol)
					{
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for(int l = 0; l < 9; ++l)
						{
							a[ll][l] -= a[icol][l]*dum;
						}
						b[ll] -= b[icol]*dum;
					}
				}
			}
		}
	}
}
	
double dateconv(double ep)
{
	if(ep > 2000000.0)
	{
		return ep - 2400000.5;
	}
	else if(ep < 4000.0) 
	{
		double nd;

		if((int)ep % 4 == 0)
		{
			nd = 366.0;
		}
		else
		{
			nd = 365.0;
		}
		return (ep-2000.0)*nd + 51544;
	}
	else
	{
		return ep;
	}
}	

int loaddata(const char *fn, Pulsar &p)
{
	FILE *in = 0;
	char str[1000];
	int l;
	char A[100], B[100], D[100];
	double ep, sr, sd;
	int n;
	const char *s;

	p.epoch = 999999.999;
	p.name[0] = 0;
	p.ref[0] = 0;
	for(int i = 5; i < 9; ++i)
	{
		p.fit[i] = 0;
	}

	l = strlen(fn);
	if(l > 128)	/* must be in-place data */
	{
		inplace = 1;
	}
	else
	{
		if(strcmp(fn, "-") == 0)
		{
			in = stdin;
		}
		else
		{
			in = fopen(fn, "r");
		}
		if(!in)
		{
			return -1;
		}
	}
	
	s = fn;
	
	for(;;)
	{
		if(inplace)
		{
			for(n = 0; s[n] != 0 && s[n] != ';'; ++n);
			if(s[n] == 0)
			{
				break;
			}
			if(n >= 999)
			{
				break;
			}
			strncpy(str, s, n);
			str[n] = 0;
			s += (n+1);
		}
		else
		{
			fgets(str, 999, in);
			if(feof(in))
			{
				break;
			}
		}
		for(int i = 0; str[i] != 0; ++i) 
		{
			if(str[i] == '#' || str[i] < ' ')
			{
				str[i] = 0;
			}
		}
		if(str[0] == 0)
		{
			continue;
		}
		if(str[0] > '9')
		{
			for(int i = 0; str[i] != 0; ++i) 
			{
				if(str[i] == '=') 
				{
					str[i] = ' ';
				}
				if(str[i] == '+' && str[i+1] == '-')
				{
					str[i] = 0;
				}
			}
			if(sscanf(str, "%s %s\n", A, B) != 2)
			{
				break;
			}
			for(int i = strlen(B) - 1; B[i] <= ' ' && i > 0; i--)
			{
				B[i] = 0;
			}
			if(strcmp(A, "RA") == 0 || strcmp(A, "ra") == 0)
			{
				p.fit[0] = 0;
				p.parm[0] = text2double(B)*PI_12;
			}
			if(strcmp(A, "Dec") == 0 || strcmp(A, "dec") == 0)
			{
				p.fit[1] = 0;
				p.parm[1] = text2double(B)*PI_180;
			}
			if(strcmp(A, "sigmu_a") == 0)
			{
				p.sigmu_a_given = atof(B);
			}
			if(strcmp(A, "sigmu_d") == 0)
			{
				p.sigmu_d_given = atof(B);
			}
			if(strcmp(A, "sigpi") == 0)
			{
				p.sigpi_given = atof(B);
			}
			if(strcmp(A, "mu_a") == 0)
			{
				p.fit[2] = 0;
				p.parm[2] = atof(B)/MAS_RAD;
			}
			if(strcmp(A, "mu_d") == 0)
			{
				p.fit[3] = 0;
				p.parm[3] = atof(B)/MAS_RAD;
			}
			if(strcmp(A, "pi") == 0 || strcmp(A, "PI") == 0 || strcmp(A, "Pi") == 0)
			{
				p.fit[4] = 0;
				p.parm[4] = atof(B)/MAS_RAD;
			}
			if(strcmp(A, "name") == 0 || strcmp(A, "Name") == 0)
			{
				strcpy(p.name, B);
			}
			if(strcmp(A, "ref") == 0 || strcmp(A, "Ref") == 0)
			{
				strcpy(p.ref, B);
			}
			if(strcmp(A, "dm") == 0 || strcmp(A, "DM") == 0)
			{
				p.dm = atof(B);
			}
			if(strcmp(A, "reverse") == 0)
			{
				p.reverse = atof(B);
				if(p.reverse > 0)
				{
					p.reverse = -1.0;
				}
				else
				{
					p.reverse = 1.0;
				}
			}
			if(strcmp(A, "break") == 0)
			{
				p.breaktime = dateconv(atof(B));
				for(int i = 5; i < 9; ++i)
				{
					p.fit[i] = 1;
				}
				fprintf(stderr, "Warning -- using obsolete feature\n");
			}
			if(strcmp(A, "sh_x") == 0)
			{
				p.fit[5] = 0;
				p.parm[5] = atof(B)/MAS_RAD;
				fprintf(stderr, "Warning -- using obsolete feature\n");
			}
			if(strcmp(A, "sh_y") == 0)
			{
				p.fit[6] = 0;
				p.parm[6] = atof(B)/MAS_RAD;
				fprintf(stderr, "Warning -- using obsolete feature\n");
			}
			if(strcmp(A, "dmu_x") == 0)
			{
				p.fit[7] = 0;
				p.parm[7] = atof(B)/MAS_RAD;
				fprintf(stderr, "Warning -- using obsolete feature\n");
			}
			if(strcmp(A, "dmu_y") == 0)
			{
				p.fit[8] = 0;
				p.parm[8] = atof(B)/MAS_RAD;
				fprintf(stderr, "Warning -- using obsolete feature\n");
			}
			if(strcmp(A, "epoch") == 0)
			{
				p.epoch = dateconv(atof(B));
			}
		}
		else
		{
			if(sscanf(str, "%lf %s %lf %s %lf",&ep,A,&sr,D,&sd)!=5)
			{
				break;
			}
			epoch[epochs] = dateconv(ep);
			a[epochs]     = text2double(A)*PI_12;
			d[epochs]     = text2double(D)*PI_180;
			siga[epochs]  = sr*PI_12/3600.0;
			sigd[epochs]  = sd*PI_180/3600.0;
			++epochs;
		}
	}

	if(strcmp(fn, "-") != 0 && inplace == 0)
	{
		fclose(in);
	}

	if(epochs == 0)
	{
		epoch[0] = 51544.0;  // year 2000.0
		a[0] = p.fit[0];
		d[0] = p.fit[1];
		siga[0] = 1.0;
		sigd[0] = 1.0;
		epochs = 1;
	}

	if(p.epoch > 999000.0)
	{
		p.epoch = epoch[0];
	}

	p.a_0 = a[0];
	p.d_0 = d[0];
	p.siga_0 = siga[0];
	p.sigd_0 = sigd[0];

	return epochs;
}

void getab(double day, double ra, double dec, double *pi_a, double *pi_d)
{
	Vector<double> alphahat, deltahat, sunhat, pulsarhat;
	Vector<double> zhat(0.0, 0.0, 1.0);

	sunhat = sunpos(day);
	pulsarhat = radec2vec(ra, dec);
	alphahat = ~(zhat^pulsarhat);
	deltahat = ~(pulsarhat^alphahat);

	*pi_d = sunhat%deltahat;
	*pi_a = sunhat%alphahat;
}

Pulsar fit(Pulsar p)
{
	Pulsar q = p;

	double fa[9], fd[9], S = 0.0, V[9], M[9][9], X[9], Minv[9][9];
	double pi_a, pi_d, W_a, W_d, pa, pd;

	for(int i = 0; i < 9; ++i)
	{
		V[i] = 0.0;
		for(int j = 0; j < 9; ++j)
		{
			M[i][j] = 0.0;
		}
	}

	for(int e = 0; e < epochs; ++e)
	{
		getab(epoch[e], p.ra_sun(epoch[e]), p.dec_sun(epoch[e]), &pi_a, &pi_d);

		pa = a[e];
		pd = d[e];

		W_a = 1.0/(siga[e]*siga[e]);
		W_d = 1.0/(sigd[e]*sigd[e]);

		fa[0] = 1.0;
		fa[1] = 0.0;
		fa[2] = (epoch[e] - p.epoch)/cos(pd)/365.242175;
		fa[3] = 0.0;
		fa[4] = pi_a/cos(pd);
		fa[5] = 0.0;
		fa[6] = 0.0;
		fa[7] = 0.0;
		fa[8] = 0.0; 

		fd[0] = 0.0;
		fd[1] = 1.0;
		fd[2] = 0.0;
		fd[3] = (epoch[e] - p.epoch)/365.242175;
		fd[4] = pi_d;
		fd[5] = 0.0;
		fd[6] = 0.0;
		fd[7] = 0.0;
		fd[8] = 0.0;
		if(epoch[e] < p.breaktime) 
		{
			fa[5] = 1.0/cos(pd);
			fd[6] = 1.0;
			fa[7] = (epoch[e] - p.epoch)/cos(pd)/365.242175;
			fd[8] = (epoch[e] - p.epoch)/365.242175;
		}

		for(int i = 0; i < 9; ++i) 
		{
			if(p.fit[i] != 1) 
			{
				pa -= fa[i]*p.parm[i];
				pd -= fd[i]*p.parm[i];
				fa[i] = fd[i] = 0.0;
			}
		}

		for(int i = 0; i < 9; ++i)
		{
			for(int j = 0; j < 9; ++j)
			{
				M[i][j] += 2.0*W_a*fa[i]*fa[j];
				M[i][j] += 2.0*W_d*fd[i]*fd[j];
			}
		}

		for(int i = 0; i < 9; ++i)
		{
			V[i] += 2.0*W_a*pa*fa[i];
			V[i] += 2.0*W_d*pd*fd[i];
		}
		S += W_a*pa*pa;
		S += W_d*pd*pd;
	}

	for(int i = 0; i < 9; ++i)
	{
		if(p.fit[i] != 1)
		{
			V[i] = p.parm[i];
			M[i][i] = 1.0;
		}
	}

	for(int i = 0; i < 9; ++i)
	{
		for(int j = 0; j < 9; ++j)
		{
			Minv[i][j] = M[i][j];
		}
	}
	for(int i = 0; i < 9; ++i)
	{
		X[i] = V[i];
	}

	gaussj9(Minv, X);

	for(int i = 0; i < 9; ++i)
	{
		if(p.fit[i] != 1) 
		{
			Minv[i][i] = 0.0;
		}
	}

	q.a_0  = X[0];
	q.d_0  = X[1];
	q.mu_a = X[2];
	q.mu_d = X[3];
	q.pi   = X[4];
	q.sh_x = X[5];
	q.sh_y = X[6];
	q.dmu_x= X[7];
	q.dmu_y= X[8];

	q.siga_0  = sqrt(Minv[0][0]);
	q.sigd_0  = sqrt(Minv[1][1]);
	q.sigmu_a = sqrt(Minv[2][2]);
	q.sigmu_d = sqrt(Minv[3][3]);
	q.sigpi   = sqrt(Minv[4][4]);
	q.sig_sh_x = sqrt(Minv[5][5]);
	q.sig_sh_y = sqrt(Minv[6][6]);
	q.sig_dmu_x= sqrt(Minv[7][7]);
	q.sig_dmu_y= sqrt(Minv[8][8]);

	for(int i = 0; i < 9; ++i)
	{
		for(int j = 0; j < 9; ++j)
		{
			q.normcov[i][j] = Minv[i][j]/sqrt(Minv[i][i]*Minv[j][j]);
		}
	}

	return q;
}

void mjd2day(double mjd, int *year, int *mon, int *day)
{
	int jd, l, n, y, m, d;

        jd = (int)(mjd+2400000);
        l = jd + 68570;
        n = 4*l/146097;
        l = l - (146097*n + 3)/4;
        y = 4000 * (l + 1) / 1461001;
        l = l - (1461*y)/4 + 31;
        m = 80*l/2447;
        d = l - (2447*m)/80;
        l = m/11;
        m = m + 2 - 12*l;
        y = 100*(n-49) + y + l;
        *year = y;
	*mon = m;
	*day = d;
}


void printday(double mjd, int opt)
{
        int y, m, d;
        char sm[12][4] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

	mjd2day(mjd, &y, &m, &d);

        printf("%d %s %02d", y, sm[m-1], d);
	if(opt == 0)
	{
		printf("\n");
	}
	else
	{
		printf(", ");
	}
}

/* Mode 0 : just print model generated from input file */

void printmodel(Pulsar *model, int argc, char **argv)
{
	int pc = 0;

	if(inplace)
	{
		model->printminimal();
	}
	else
	{
		model->print();
		for(int i = 2; i < argc; ++i)
		{
			double y;
			
			if(strcmp(argv[i], "-c") == 0 && pc == 0) 
			{
				model->printnormcov();
				pc = 1;
			}
			y = atof(argv[i]);
			if(y > 0.0)
			{
				char ras[100], decs[100];

				if(y < 1950.0 || y > 2050.0)
				{
					continue;
				}
				double2text(model->ra_earth(y)/PI_12, ras);
				double2text(model->dec_earth(y)/PI_180, decs);
				printf("%8.3f %s %s\n", y, ras, decs);
			}
		}
	}
}

/* Mode 5 : write psr file */

void makepsrfile(Pulsar *model, int argc, char **argv)
{
	model->saveaspsr();
}

/* Mode 1 : Predictor mode, earth coords */

void predictearth(Pulsar *model, int argc, char **argv)
{
	double start, stop, f, inc, sigx=0.0, sigy=0.0;

	if(argc != 4 && argc != 6 && argc != 8)
	{
		fprintf(stderr, "Wrong number of arguments, buddy\n");
		fprintf(stderr, "Run with no arguments for help\n");

		return;
	}

	start = atof(argv[3]);
	if(start < 4000.0)
	{
		f = 365.242175;
	}
	else
	{
		f = 1.0;
	}
	start = dateconv(start);
	if(argc == 4)
	{
		inc = 2000.0*f;
		stop = start;
	}
	else
	{
		stop = dateconv(atof(argv[4]));
		inc = atof(argv[5])*f;
	}
	if(inc <= 0.0)
	{
		fprintf(stderr, "Error -z : Need positive increment\n");
		
		return;
	}
	if(argc == 8)
	{
		sigx = atof(argv[6]);
		sigy = atof(argv[7]);
	}

	for(double mjd = start; mjd <= stop; mjd+=inc)
	{
		char ras[100], decs[100];

		double2text(model->ra_earth(mjd)/PI_12, ras, ':');
		double2text(model->dec_earth(mjd)/PI_180,decs, ':');
		if(argc == 8)
		{
			printf("%9.3f  %s %f  %s %f\n", mjd, ras, sigx, decs, sigy);
		}
		else
		{
			printf("%9.3f  %s  %s\n", mjd, ras, decs);
		}
	}
}

/* Mode 2 : Predictor mode, sun coords */

void predictsun(Pulsar *model, int argc, char **argv)
{
	double start, stop, inc, f;

	if(argc != 4 && argc != 6)
	{
		printf("Wrong number of arguments, buddy\n");
		printf("Run with no arguments for help\n");
		
		return;
	}

	start = atof(argv[3]);
	if(start < 4000.0)
	{
		f = 365.242175;
	}
	else
	{
		f = 1.0;
	}
	start = dateconv(start);
	if(argc == 4)
	{
		inc = 2000.0*f;
		stop = start;
	}
	else
	{
		stop = dateconv(atof(argv[4]));
		inc = atof(argv[5])*f;
	}
	if(inc <= 0.0)
	{
		fprintf(stderr, "Error -s : Need positive increment\n");
		
		return;
	}

	for(double mjd = start; mjd <= stop; mjd+=inc)
	{
		char ras[100], decs[100];
		
		double2text(model->ra_sun(mjd)/PI_12, ras, ' ');
		double2text(model->dec_sun(mjd)/PI_180,decs,' ');
		printf("%9.3f  %s  %s\n", mjd, ras, decs);
	}
}

/* Mode 3 : Find parallax extrema */

void findextrema(Pulsar *model, double center)
{
	double maxa = -1.0, mina = 1.0;
	double maxd = -1.0, mind = 1.0;
	double maxt = -1.0, mint = 1.0;
	double maxt2 = -1.0, mint2 = 1.0;
	double tmaxa = 0.0, tmina = 0.0;
	double tmaxd = 0.0, tmind = 0.0;
	double tmaxt = 0.0, tmint = 0.0;
	double tmaxt2 = 0.0, tmint2 = 0.0;
	double ra, dec;

	ra = model->ra_sun(center);
	dec = model->dec_sun(center);
	if(center < 0.0)
	{
		center = model->epoch;
	}

	for(double mjd = center-182.5; mjd <= center+182.5; mjd+=0.25)
	{
		double pi_a, pi_d, pi_t;

		getab(mjd, ra, dec, &pi_a, &pi_d);
		pi_t = sqrt(pi_a*pi_a + pi_d*pi_d);
		if(pi_a < mina) 
		{ 
			mina = pi_a; 
			tmina = mjd;
		}
		if(pi_a > maxa)
		{
			maxa = pi_a;
			tmaxa = mjd;
		}
		if(pi_d < mind)
		{
			mind = pi_d;
			tmind = mjd;
		}
		if(pi_d > maxd)
		{
			maxd = pi_d;
			tmaxd = mjd;
		}
		if(mjd < center)
		{
			if(pi_t < mint)
			{
				mint = pi_t;
				tmint = mjd;
			}
			if(pi_t > maxt)
			{
				maxt = pi_t;
				tmaxt = mjd;
			}
		}
		else
		{
		 	if(pi_t < mint2)
			{
				mint2 = pi_t;
				tmint2 = mjd;
			}
		 	if(pi_t > maxt2)
			{
				maxt2 = pi_t;
				tmaxt2 = mjd;
			}
		}
	}
	printf("Maximum pi in RA is %5.3f occuring at %9.3f = ", maxa*model->pi*MAS_RAD, tmaxa);
	printday(tmaxa);
	printf("Minimum pi in RA is %5.3f occuring at %9.3f = ", mina*model->pi*MAS_RAD, tmina);
	printday(tmina);
	printf("Maximum pi in Dec is %5.3f occuring at %9.3f = ", maxd*model->pi*MAS_RAD, tmaxd);
	printday(tmaxd);
	printf("Minimum pi in Dec is %5.3f occuring at %9.3f = ", mind*model->pi*MAS_RAD, tmind);
	printday(tmind);
	printf("Maximum total pi is %5.3f occuring at %9.3f and %9.3f = ", maxt*model->pi*MAS_RAD, tmaxt, tmaxt2);
	printday(tmaxt, 1);
	printday(tmaxt2);
	printf("Minimum total pi is %5.3f occuring at %9.3f and %9.3f = ", mint*model->pi*MAS_RAD, tmint, tmint2);
	printday(tmint, 1);
	printday(tmint2);
}
			
/* Mode 4 : produce output plots */

void makeplots(Pulsar *model, int dropm, int dropp, int drops)
{
	FILE *outx;
	double x0 = 0.0, y0 = 0.0;
	double ee, extra = 0.25;

	ee = epoch[epochs-1] - epoch[0];
	if(ee + 2*extra < 1.2)
	{
		extra = 0.5*(1.4-ee);
	}

	outx = fopen("pmpar_t", "w");
	for(double t = epoch[0] - extra; t < epoch[epochs-1] + extra; t += 1.0)
	{
		double x, y;

		model->maspos(t, dropm, dropp, drops, &x, &y);
		if(t > model->breaktime && t-1.0 <= model->breaktime)
		{
			fprintf(outx, "\n");
		}
		fprintf(outx, "%f %f %f\n", t, x-x0, y-y0);
	}
	fclose(outx);

	outx = fopen("pmpar_e", "w");
	for(int e = 0; e < epochs; ++e)
	{
		double x, y, xm, ym;

		model->maspos(epoch[e], dropm, dropp, drops, &xm, &ym);
		model->maspoint(e, dropm, dropp, drops, &x, &y);
		fprintf(outx, "%f %f %f %f %f %f %f\n", 
			epoch[e], x-x0, siga[e]*MAS_RAD*cos(d[e]), 
			y-y0, sigd[e]*MAS_RAD, xm-x0, ym-y0);
	}
	fclose(outx);
}


int main(int argc, char **argv)
{
	Pulsar model;
	int mode = 0;
	double center = -1;
	int dropm = 0, dropp = 0, drops = 0;

	if(argc < 2)
	{
		printf("pmpar version %s  %s\n", VERSION, PACKAGE_BUGREPORT);
		printf("A parallax and proper motion fitting machine\n\n");
		printf("Usage : pmpar <input file> [epoch1] [epoch2] ...\n\n");
		printf("  Where <input file> has the following format:\n");
		printf("          [format to de described later]\n");
		printf("  And epoch1, epoch2, ... are years to compute position for.\n\n");
		printf("        pmpar <input file> -{z|s} <start> [<stop> <inc>]\n\n");
		printf("  is predictor mode.  -z is for earth coords, -s for barycenter coords\n");
		printf("  start, stop, inc in years, for ephemeris to be made\n\n");
		printf("        pmpar <input file> -m [epoch]\n\n");
		printf("  determines min and max parallax, and when.  Epoch defaults to epoch in file\n\n");
		printf("        pmpar <input file> -o[m][p][s]\n\n");
		printf("  generates output files pmpar_t and pmpar_e\n\n");
		printf("  Please read pmpar.doc for a better description\n");

		return 0;
	}

// Syncronize model with data

	if(loaddata(argv[1], model) < 0)
	{
		return -1;
	}

	for(int iter = 0; iter < 3; ++iter)
	{
		model = fit(model);
	}

	if(argc > 2) 
	{
		if(strcmp(argv[2], "-z") == 0)
		{
			mode = 1;
		}
		if(strcmp(argv[2], "-s") == 0)
		{
			mode = 2;
		}
		if(strcmp(argv[2], "-p") == 0)
		{
			mode = 5;
		}
		if(strcmp(argv[2], "-m") == 0)
		{
			mode = 3;
			if(argc > 3)
			{
				center = atof(argv[3]);
			}
		}
		if(argv[2][0] == '-' && argv[2][1] == 'o')
		{
			mode = 4;
			for(unsigned int r = 2; r <= 4 ; r++)
			{
				if(strlen(argv[2]) > r)
				{
					if(argv[2][r] == 'm') dropm = 1;
					if(argv[2][r] == 'p') dropp = 1;
					if(argv[2][r] == 's') drops = 1;
				}
			}
		}
	}

	if(!inplace)
	{
		makeplots(&model, dropm, dropp, drops);
	}

	switch(mode)
	{
	case 0:
		printmodel(&model, argc, argv);
		break;
	case 1:
		predictearth(&model, argc, argv);
		break;
	case 2:
		predictsun(&model, argc, argv);
		break;
	case 3:
		findextrema(&model, center);
		break;
	case 5:
		makepsrfile(&model, argc, argv);
		break;
	}

	return 0;
}
