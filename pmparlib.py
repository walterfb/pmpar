#!/usr/bin/env python

from math import *
from libplot import Plotter
from string import strip, split, find
import biggles
import os
import sys
from pdf import PDF

def get1sigma(file):
	p = PDF(file)
	pk = p.peak()
	sig1 = p.shortestinterval(0.68268949)
	return pk, sig1[0], sig1[1]


class pmparpsr:
	def __init__(self, filename):
		if os.access(filename, os.F_OK) == 0:
			self.name = "None"
			self.filename = filename
			return
		output = os.popen("pmpar %s" % filename, "r").readlines();
		for l in output:
			s = split(strip(l));
			if len(s) < 3:
				continue
			if s[0] == 'Name':
				self.name = s[2]
			if s[0] == 'l':
				self.l = float(s[2])*3.1415926535/180.0
			if s[0] == 'b':
				self.b = float(s[2])*3.1415926535/180.0
			if len(s) < 5:
				continue
			if s[0] == 'mu_a':
				self.mu_a = float(s[2])
				self.dmu_a = float(s[4])
			if s[0] == 'mu_d':
				self.mu_d = float(s[2])
				self.dmu_d = float(s[4])
			if s[0] == 'pi':
				self.pi = float(s[2])
				self.dpi = float(s[4])
			if s[0] == 'ne':
				self.ne = float(s[2])
				self.dne = float(s[4])
#			if s[0] == 'dist':
#				self.dist = float(s[2])
#				self.distmin = self.dist - float(s[4])
#				if self.distmin < 0.0:
#					self.distmin = 0.0
#				self.distmax = self.dist + float(s[6])
		(self.dist, self.distmin, self.distmax) = get1sigma(split(filename,'.')[0]+'.dist')
		

	def prt(self):
		print self.name, ":", self.mu_a, "+-", self.dmu_a, ",", self.mu_d, "+-",self.dmu_d, ",", self.pi,"+-", self.dpi

def multipmpar(filelist):
	psrs = []
	for f in filelist:
		psrs.append(pmparpsr(f))
	return psrs

if __name__ == '__main__':
	psrs = multipmpar(sys.argv[1:])
	for p in psrs:
		p.prt()
