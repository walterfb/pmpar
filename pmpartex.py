#!/usr/bin/env python

import sys
import os
from string import strip, split

def texify(filename):
	command = "pmpar " + filename
	name = ''
	ref = ''
	ra = ''
	dec = ''
	mua = ''
	mud = ''
	pi = ''
	output = os.popen(command, "r").readlines()
	for o in output:
		l = strip(o)
		if l[0:4] == 'Name':
			name = strip(l[7:])
		if l[0:4] == 'Ref.':
			ref = strip(l[7:])
		if l[0:2] == 'RA':
			p = split(split(l[7:])[0], ':')
			d = split(p[2], '.')
			ra = "$"+p[0]+"^h"+p[1]+"^m"+d[0]+"^s."+d[1][0:4]+"$"
		if l[0:3] == 'Dec':
			p = split(split(l[7:])[0], ':')
			d = split(p[2], '.')
			dec = "$"+p[0]+"^{\\circ}"+p[1]+"'"+d[0]+"''."+d[1][0:4]+"$"
		if l[0:4] == 'mu_a':
			mu = split(l[7:])
			mua = ("$%4.2f"%float(mu[0]))+r"\pm" \
				+("%4.2f$"%float(mu[2]))
		if l[0:4] == 'mu_d':
			mu = split(l[7:])
			mud = ("$%4.2f"%float(mu[0]))+r"\pm" \
				+("%4.2f$"%float(mu[2]))
		if l[0:3] == 'pi ':
			p = split(l[7:])
			pi = ("$%4.2f"%float(p[0]))+r"\pm" \
				+("%4.2f$"%float(p[2]))
	print name,"&",ra,"&",dec,"&",mua,"&",mud,"&",pi,r"\\"

def texifyall(filenames):
	for f in filenames:
		texify(f)


filenames = []

for a in sys.argv[1:]:
	if a[0] != '-':
		filenames.append(a)
texifyall(filenames)
