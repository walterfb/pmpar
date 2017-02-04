#!/usr/bin/env python

from math import *
from libplot import Plotter
from string import strip, split, find
import biggles
import os
import sys
from pmparlib import *

def makeneplot(filelist):
	psrs = multipmpar(filelist)
	X = []
	Y = []
	p = biggles.Plot()
	p.xrange = -3.8,3.8
	p.yrange = -3.8,3.8
	p.add(biggles.Line((-3.2, 0.0), (3.2, 0.0)))
	p.add(biggles.Line((0.0 ,-3.2), (0.0, 3.2)))
	p.add(biggles.Circle(0.0, 0.0, 1.0, linewidth = 0.5, linetype='dotted'))
	p.add(biggles.Circle(0.0, 0.0, 2.0, linewidth = 0.5, linetype='dotted'))
	p.add(biggles.Circle(0.0, 0.0, 3.0, linewidth = 0.5, linetype='dotted'))
	p.add(biggles.DataLabel(0.0, -1.0, "1 kpc", halign = "left", valign = "top", fontsize = 2))
	p.add(biggles.DataLabel(0.0, -2.0, "2 kpc", halign = "left", valign = "top", fontsize = 2))
	p.add(biggles.DataLabel(0.0, -3.0, "3 kpc", halign = "left", valign = "top", fontsize = 2))
	p.add(biggles.DataLabel(3.4, 0.0, r"$\ell = 0$", haligh = "left", valighn = "center", fontsize = 2))
	p.add(biggles.DataLabel(0.0, 3.2, r"$\ell = 90$", haligh = "center", valign = "bottom", fontsize = 2))
	lva = "center"
	la = None
	deltay = 0.0
	delta = 0.1
	for psr in psrs:
		if psr.name == "None":
			if psr.filename == '+':
				lva = "bottom"
			if psr.filename == '-':
				lva = "top"
			if psr.filename == 'l':
				la = "left"
			if psr.filename == 'r':
				la = "right"
			if psr.filename == 'L':
				la = "left"
				delta = delta + 0.12
			if psr.filename == 'R':
				la = "right"
				delta = delta + 0.12
			if psr.filename == '^':
				deltay = deltay+0.066
			if psr.filename == '_':
				deltay = deltay-0.066
			continue
		x = cos(psr.l) #* cos(psr.b)
		y = sin(psr.l) #* cos(psr.b)
		p.add(biggles.Line( (x*psr.distmin, y*psr.distmin), (x*psr.distmax, y*psr.distmax) ))
		X.append(x*psr.dist)
		Y.append(y*psr.dist)
		print psr.name, psr.dist
		if la == None:
			if x > 0:
				la = "left"
			else:
				la = "right"
		if la == 'left':
			lx = x*psr.dist + delta
		else:
			lx = x*psr.dist - delta
		ly = y*psr.dist + deltay
#		p.add(biggles.DataLabel(lx, ly, psr.name, halign = la, valign = "bottom", fontsize = 1.5))
#		p.add(biggles.DataLabel(lx, ly, "$n_e$=%6.4f" % psr.ne, halign = la, valign = "top", fontsize = 1.5))
		p.add(biggles.DataLabel(lx, ly, psr.name+" $n_e$=%6.4f" % psr.ne, halign = la, valign = lva, fontsize = 1.4))
		lva = "center"
		delta = 0.1
		deltay = 0.0
		la = None

	rad = [0.05]*len(X)

#	p.add(biggles.Points(X, Y, symboltype="circle" , size=1.25) )
#	p.add(biggles.Points(X, Y, symboltype="circle" , size=1.5) )
	p.add(biggles.Circles(X, Y, rad, width=1.6))

	return p

p = makeneplot(sys.argv[1:])
p.show()
p.save_as_eps("neplot.ps")
