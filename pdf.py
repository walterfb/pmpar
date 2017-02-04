#!/usr/bin/env python

from string import split, strip

# the data are stored a bit differently in the python implementation than in
# the C implementation to make plotting easier.
# This is basically a direct translation of a subset of pdf.c into python.

class PDF:
	def load(self, filename):
		self.v = []
		self.pdf = []
		lines = open(filename).readlines()
		for l in lines[1:]:
			data = split(strip(l))
			self.v.append(float(data[0]))
			self.pdf.append(float(data[1]))
		self.dv = (self.v[-1] - self.v[0])/(len(self.v) - 1.0) 

	def save(self, filename):
		f = open(filename, 'w')
		f.write('# %f %d\n' % (self.v[0], len(self.pdf)))
		for i in range(0, len(self.pdf)):
			f.write('%f %f\n' % (self.v[i], self.pdf[i]))
		f.close()
			
	def __init__(self, filename=None):
		self.v = []
		self.pdf = []
		self.dv = 0
		if filename != None:
			self.load(filename)

	def accumulate(self):
		pdf = PDF()
		pdf.dv = self.dv;
		pdf.pdf.append(self.pdf[0]*self.dv)
		pdf.v.append(self.v[0])
		for i in range(1, len(self.v)):
			pdf.v.append(self.v[i])
			pdf.pdf.append(pdf.pdf[i-1] + self.dv*self.pdf[i])
		return pdf

	def peak(self, v1=None, v2=None):
		p = -9999.999
		q = 0.0
		if v1 == None:
			v1 = self.v[0]
		if v2 == None:
			v2 = self.v[-1]
		for i in range(0, len(self.v)):
			if self.v[i] >= v1 and self.v[i] <= v2:
				if self.pdf[i] > p:
					p = self.pdf[i]
					q = self.v[i]
		return q;

	def totalarea(self):
		sum = 0.0
		for p in self.pdf:
			sum = sum + p
		return sum*self.dv

	def shortestinterval(self, area, v1=None, v2=None):
		goal = area*self.totalarea()/self.dv
		i1 = -1
		i2 = -1
		if v1 == None:
			v1 = self.v[0]
		if v2 == None:
			v2 = self.v[-1]
		for i in range(0, len(self.v)):
			if i1 < 0 and self.v[i] >= v1:
				i1 = i
			if i2 < 0 and self.v[i] > v2:
				i2 = i-1
		if i2 < 0:
			i2 = len(self.v) - 1
		if i1 < 0 or i2 < 0:
			return (0.0, 0.0)
		j1 = j2 = i1
		b1 = b2 = -1
		bsize = len(self.v)*2
		sum = self.pdf[i1]
		while(j2 <= i2):
			if sum < goal or j1 >= j2:
				j2 = j2+1
				if j2 > i2:
					break
				sum = sum + self.pdf[j2]
				if sum >= goal and j2-j1 < bsize:
					b1, b2 = j1, j2
					bsize = b2-b1
			else:
				sum = sum - self.pdf[j1]
				if sum > goal and j2-j1 < bsize:
					b1, b2 = j1, j2
					bsize = b2-b1
				j1 = j1+1
		if b1 < 0:
			return (0.0, 0.0)
		return (self.v[b1], self.v[b2])
		
	def scale(self, factor):
		for i in range(0, len(self.pdf)):
			self.pdf[i] = self.pdf[i]*factor

	def norm(self):
		sum = 0.0
		for p in self.pdf:
			sum = sum + p
		if sum != 0.0:
			self.scale(1.0/(sum*self.dv))

	def interp(self, v):
		index = int((v - self.v[0])/self.dv)
		if (index < 0 or index >= len(self.pdf)): 
			a = 0.0
		else:
			a = self.pdf[index]

		index = index+1
		if (index < 0 or index >= len(self.pdf)):
			b = 0.0
		else:
			b = self.pdf[index]

		f = float(index) - (v - self.v[0])/self.dv;

		return f*a + (1.0-f)*b
