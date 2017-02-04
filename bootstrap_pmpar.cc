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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "config.h"

#define MAXEPOCHS	100

typedef struct
{
	double mjd, ra, dec, dra, ddec;
} Datum;

void writeDatum(FILE *out, const Datum *d)
{
	fprintf(out, "%12.6f %14.7f %9.7f %14.7f %9.7f\n", d->mjd, d->ra, d->dra, d->dec, d->ddec);
}

/* only loads matching name, ifnum */
int parseDatum(Datum *d, char *name, const char *str, int ifnum)
{
	int i, n, epo;
	char psr[50], expname[50];
	double dra1, dra2, ddec1, ddec2;
	double F;

	n = sscanf(str, "%s %s %lf %d %d %lf %lf %lf %lf %lf %lf",
		psr, expname, &d->mjd, &epo, &i, 
		&d->ra,  &dra1,  &dra2,
		&d->dec, &ddec1, &ddec2);
		
	if(n != 11) 
	{
		return 0;
	}

	if(ifnum != i && ifnum != -1)
	{
		return 0;
	}

	if(name[0] == 0)
	{
		strcpy(name, psr);
	}
	else if(strcmp(name, psr) != 0)
	{
		return 0;
	}

	/* FIXME -- verify assignment of errors is correct here */

	F = 15.0*cos(d->dec*M_PI/(3600.0*180.0));
	
	d->dra  = dra1/1000.0;
	d->dra /= F;
	d->ddec = ddec1/1000.0;

	return 1;
}

int loaddata(Datum *d, char *name, const char *filename, int ifnum)
{
	FILE *in;
	char line[257];
	int n;

	in = fopen(filename, "r");
	if(!in)
	{
		fprintf(stderr, "Cannot open %s for read\n", filename);
		return -1;
	}

	n = 0;
	for(;;)
	{
		fgets(line, 256, in);
		if(feof(in)) break;
		n += parseDatum(d+n, name, line, ifnum);
	}
	fclose(in);

	if(n < 3)
	{
		fprintf(stderr, "Not enough epochs found for solution\n");
	}

	return n;
}

void bootstrap_pmpar(FILE *out, const Datum *data, int epochs, int nonrandom)
{
	FILE *p;
	char str[1000];
	char command[32768];
	const Datum *d;
	int i, r;
	char epochstr[4*MAXEPOCHS+1];

	strcpy(command, "pmpar \"epoch=2002;");

	epochstr[0] = 0;

	for(i = 0; i < epochs; ++i)
	{
		if(nonrandom)
		{
			r = i;
		}
		else
		{
			r = (lrand48() / 17) % epochs;
		}

		sprintf(str, ",%d", r);
		strcat(epochstr, str);

		d = data+r;
	
		sprintf(str, "%12.6f %14.7f %9.7f %14.7f %9.7f;", 
			d->mjd, d->ra, d->dra, d->dec, d->ddec);

		strcat(command, str);
	}

	epochstr[0] = ' ';

	strcat(command, "\"");

	p = popen(command, "r");

	if(!p)
	{
		fprintf(stderr, "cannot open pipe to pmpar\n");

		exit(0);
	}
	
	fgets(str, 999, p);
	for(i = 0; str[i]; i++)
	{
		if(str[i] < ' ')
		{
			str[i] = ' ';
		}
	}
	strcat(str, epochstr);
	strcat(str, "\n");
	fclose(p);
	fprintf(out, "%s", str);
}

int usage(const char *program)
{
	fprintf(stderr, "bootstrap_pmpar version %s -- %s\n\n", VERSION, PACKAGE_BUGREPORT);
	fprintf(stderr, "Usage : %s <input file> <ifnum> [<num iterations>]\n\n", program);
	fprintf(stderr, "  <input file> is of .dat format\n");
	fprintf(stderr, "  <ifnum> is integer between 1 and 16, or \"all\" for all IFs\n");
	fprintf(stderr, "  <num interations> is integer > 0 [optional -- default = 1000]\n\n");
	fprintf(stderr, "Output is sent to stdout.  The 11 columns of output are:\n");
	fprintf(stderr, "  1.  RA (in seconds of time)\n");
	fprintf(stderr, "  2.  sigma RA (in seconds of time)\n");
	fprintf(stderr, "  3.  Dec (in arcsec)\n");
	fprintf(stderr, "  4.  sigma Dec (in arcsec)\n");
	fprintf(stderr, "  5.  mu_a (in mas/yr)\n");
	fprintf(stderr, "  6.  sigma mu_a (in mas/yr)\n");
	fprintf(stderr, "  7.  mu_d (in mas/yr)\n");
	fprintf(stderr, "  8.  sigma mu_d (in mas/yr)\n");
	fprintf(stderr, "  9.  pi (in mas)\n");
	fprintf(stderr, "  10. sigma pi (in mas)\n");
	fprintf(stderr, "  11. epoch mask -- list of epoch numbers used (0-based)\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note that pmpar version 1.6 or greater must be in your path\n");

	return 0;
}

int main(int argc, char **argv)
{
	Datum data[MAXEPOCHS];
	FILE *out;
	char psr[100];
	int ifnum = 1;
	int epochs;
	int i;
	int nboot = 1000;
	
	srand48(time(0));

	if(argc < 3)
	{
		return usage(argv[0]);
	}

	if(strcmp(argv[2], "all") == 0)
	{
		ifnum = -1;
		fprintf(stderr, "All IFs being used together\n");
	}
	else
	{
		ifnum = atoi(argv[2]);

		if(ifnum <= 0 || ifnum > 8)
		{
			return usage(argv[0]);
		}
	}

	if(argc > 3)
	{
		nboot = atoi(argv[3]);
		if(nboot <= 0)
		{
			return usage(argv[0]);
		}
	}
	
	epochs = loaddata(data, psr, argv[1], ifnum);

	if(epochs < 3)
	{
		fprintf(stderr, "Sorry, data not loaded\n");

		return 0;
	}

	fprintf(stderr, "%d epochs of data being used\n", epochs);

	out = stdout;

	bootstrap_pmpar(out, data, epochs, 1);
		
	for(i = 0; i < nboot; i++)
	{
		bootstrap_pmpar(out, data, epochs, 0);
	}

	return 0;
}
