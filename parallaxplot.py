#!/usr/bin/env python

from math import *
from string import strip, split, find
import biggles
import os
import sys

def autosetcustomticks(p):
	Dx = p.xrange[1] - p.xrange[0]
	Dy = p.yrange[1] - p.yrange[0]
	if Dx > Dy:
		D = Dx
	else:
		D = Dy
	fac=4
	dmit = 0.25
	if D > 9.0:
		fac=4
		dmit = 0.5
	if D > 21.0:
		fac=5
		dmit = 1.0
	if D > 35.0:
		fac=4
		dmit = 2.5
	if D > 90.0:
		fac=4
		dmit = 5.0
	if D > 210.0:
		fac=5
		dmit = 10.0
	if D > 350.0:
		fac=4
		dmit = 25.0
	if D > 900.0:
		fac=4
		dmit = 50.0
	if D > 2100.0:
		fac=5
		dmit = 100.0
	p.x1.major_ticks = p.x2.major_ticks = []
	p.x1.minor_ticks = p.x2.minor_ticks = []
	p.x2.major_ticks = p.x2.major_ticks = []
	p.x2.minor_ticks = p.x2.minor_ticks = []
	p.x1.major_ticklabels = []
	p.y.major_ticks = p.y2.major_ticks = []
	p.y.minor_ticks = p.y2.minor_ticks = []
	min = int(p.xrange[0]/dmit)
	max = int(p.xrange[1]/dmit)
	if min*dmit < p.xrange[0]:
		min = min+1
	if max*dmit > p.xrange[1]:
		max = max-1
	for i in range(min, max+1):
		if i%fac == 0:
			p.x1.major_ticks.append(i*dmit)
			p.x2.major_ticks.append(i*dmit)
			p.x1.major_ticklabels.append("%d"%(-i*dmit))
		else:
			p.x.minor_ticks.append(i*dmit)
	min = int(p.yrange[0]/dmit)
	max = int(p.yrange[1]/dmit)
	if min*dmit < p.yrange[0]:
		min = min+1
	if max*dmit > p.yrange[1]:
		max = max-1
	for i in range(min, max+1):
		if i%fac == 0:
			p.y.major_ticks.append(i*dmit)
		else:
			p.y.minor_ticks.append(i*dmit)

def makeparallaxplot(argv):
	name = 'Pulsar'

	if find(sys.argv[-1], "C") > 0:
		docolor = 1
		curvecolor = 'green'
		ellipsecolor = 'black'
		stringcolor = 'black'
	else:
		docolor = 0
		curvecolor = 'black'
		ellipsecolor = 'black'
		stringcolor = 'black'
	if find(sys.argv[-1], "a") > 0:
		fsize = 3.4
		fspace = 0.06
	else:
		fsize = 2.2
		fspace = 0.03

	if find(sys.argv[-1], "l") > 0:
		dolabel = 1
	else:
		dolabel = 0
	if find(sys.argv[-1], "b") > 0:
		dobars = 1
	else:
		dobars = 0

	if find(sys.argv[-1], "d") > 0:
		dodots = 1
	else:
		dodots = 0

	if find(sys.argv[-1], "m") > 0:
		mtrue = 1
	else:
		mtrue = 0

	cmd = "pmpar %s" % argv[1]
	print 'executing ', cmd
	output = os.popen(cmd, "r").readlines();
	for l in output:
		s = split(strip(l));
		if len(s) < 3:
			continue
		if s[0] == 'Name':
			name = s[2]
		if len(s) < 5:
			continue
		if s[0] == 'mu_a':
			mu_a = float(s[2])
			dmu_a = float(s[4])
		if s[0] == 'mu_d':
			mu_d = float(s[2])
			dmu_d = float(s[4])
		if s[0] == 'pi':
			pi = float(s[2])
			dpi = float(s[4])
		
	mu_a_str = "$\\mu_{\\alpha}$=%4.2f$\\pm$%4.2f mas/yr" % (mu_a, dmu_a)
	mu_d_str = "$\\mu_{\\delta}$=%4.2f$\\pm$%4.2f mas/yr" % (mu_d, dmu_d)
	pi_str   = "$\\pi$=%4.2f$\\pm$%4.2f mas" % (pi  , dpi  )

	if(len(argv) > 2):
		os.system("pmpar %s %s" % (argv[1], argv[2]));
	else:
		os.system("pmpar %s -o" % argv[1]);

	txdata = []
	tydata = []
	t0 = -100.0
	d1 = open("pmpar_t", "r").readlines()
	for d in d1:
		s = split(strip(d))
		if(len(s) > 2):
			if t0 < 0.0:
				t0 = float(s[0])
			if (float(s[0]) > t0 + 1.0 and mtrue == 1):
				break;
			txdata.append(-float(s[1]))
			tydata.append(float(s[2]))

	if mtrue == 1:
		txdata.append(txdata[0])
		tydata.append(tydata[0])

	m1 = max(txdata)
	m2 = max(tydata)
	n1 = min(txdata)
	n2 = min(tydata)

	p = biggles.FramedPlot()
	
	exdata = []
	ex1data = []
	ex2data = []
	eydata = []
	ey1data = []
	ey2data = []
	erxdata = []
	erydata = []
	mxdata = []
	mydata = []
	tdata = []
	d2 = open("pmpar_e", "r").readlines()
	for d in d2:
		s = split(strip(d))
		if len(s) > 4:
			tdata.append(float(s[0]))
			exdata.append(-float(s[1]))
			eydata.append(float(s[3]))
			erxdata.append(float(s[2]))
			erydata.append(float(s[4]))
			ex1data.append(-float(s[1])-float(s[2]))
			ex2data.append(-float(s[1])+float(s[2]))
			ey1data.append(float(s[3])-float(s[4]))
			ey2data.append(float(s[3])+float(s[4]))
			mxdata.append(-float(s[5]))
			mydata.append(float(s[6]))
	
	m = min(min(ex1data), -2.0)
	if m < n1 : n1 = m
	m = max(max(ex2data), 2.0)
	if m > m1 : m1 = m
	m = min(min(ey1data), -2.0)
	if m < n2 : n2 = m
	m = max(max(ey2data), 2.0)
	if m > m2 : m2 = m

	m = m*1.1

	r1 = m1 - n1
	r2 = m2 - n2
	if r1 > r2:
		m2 = m2 + 0.5*(r1-r2)
		n2 = n2 - 0.5*(r1-r2)
	else:
		m1 = m1 + 0.5*(r2-r1)
		n1 = n1 - 0.5*(r2-r1)
	r = m1 - n1
	
	p.xrange = n1-0.05*r, m1+0.05*r
	p.yrange = n2-0.05*r, m2+0.05*r
	p.title = "Motion of %s" % name
	p.xlabel = r"$\Delta\alpha$ [mas]"
	p.ylabel = r"$\Delta\delta$ [mas]"
	if dodots == 1:
		p.add(biggles.Points(mxdata, mydata, \
			symboltype="filled circle", symbolsize = 1.0))
	if dobars == 1:
		for i in range(0, len(mxdata)):
			p.add(biggles.Line((mxdata[i], mydata[i]), \
				(exdata[i], eydata[i]), linewidth=0.5))
	p.add(biggles.Curve(txdata, tydata, linewidth=0.8, color=curvecolor))
	p.add(biggles.Ellipses(exdata, eydata, erxdata, erydata, linewidth=0.5,color=ellipsecolor))
#	p.add(biggles.XErrorBars(eydata, ex1data, ex2data))
#	p.add(biggles.YErrorBars(exdata, ey1data, ey2data))
	if dolabel == 1:
		for i in range(0, len(tdata)):
			t = tdata[i]
			if mtrue == 0:
				hal = "left"
				val = "center"
				f = 2.5
				if(t < 1998.335):
					f = -2.5
					hal = "right"
			else:
				if exdata[i] < 0:
					hal = "left"
					f = 1.3
				else:
					hal = "right"
					f = -1.3
				if eydata[i] < 0:
					val = "bottom"
				else:
					val = "top"
				
			x = exdata[i] + f*erxdata[i]
			y = eydata[i]
			p.add(biggles.DataLabel(x, y, "%8.3f" % t, \
				valign = val, halign = hal, fontsize = fsize))
	p.comment = name
	
	if find(sys.argv[-1], "S") > 0:
		p.add(biggles.PlotLabel(0.05, 1.0-fspace, mu_a_str, fontsize = fsize, \
			valign = 'top', halign = 'left', color=stringcolor))
		p.add(biggles.PlotLabel(0.05, 1.0-2.0*fspace, mu_d_str, fontsize = fsize, \
			valign = 'top', halign = 'left', color=stringcolor))
		p.add(biggles.PlotLabel(0.05, 1.0-3.0*fspace, pi_str, fontsize = fsize, \
			valign = 'top', halign = 'left', color=stringcolor))

	return p

def makeparallaxplotxy(argv):
	name = 'Pulsar'

	if find(sys.argv[-1], "q") > 0:
		showsigpi = 1
	else:
		showsigpi = 0
	output = os.popen("pmpar %s" % argv[1], "r").readlines();
	for l in output:
		s = split(strip(l));
		if len(s) < 3:
			continue
		if s[0] == 'Name':
			name = s[2]
		if len(s) < 5:
			continue
		if s[0] == 'mu_a':
			mu_a = float(s[2])
			dmu_a = float(s[4])
		if s[0] == 'mu_d':
			mu_d = float(s[2])
			dmu_d = float(s[4])
		if s[0] == 'pi':
			pi = float(s[2])
			dpi = float(s[4])
		
	mu_a_str = "$\\mu_{\\alpha}$=%4.2f$\\pm$%4.2f mas/yr" % (mu_a, dmu_a)
	mu_d_str = "$\\mu_{\\delta}$=%4.2f$\\pm$%4.2f mas/yr" % (mu_d, dmu_d)
	pi_str   = "$\\pi$=%4.2f$\\pm$%4.2f mas" % (pi  , dpi  )

	if(len(argv) > 2):
		os.system("pmpar %s %s" % (argv[1], argv[2]));
	else:
		os.system("pmpar %s -o" % argv[1]);

	ttdata = []
	txdata = []
	txdatam = []
	txdatap = []
	tydata = []
	tydatam = []
	tydatap = []
	fm = (pi-dpi)/pi
	fp = (pi+dpi)/pi
	d1 = open("pmpar_t", "r").readlines()
	for d in d1:
		s = split(strip(d))
		if(len(s) > 2):
			txdata.append(float(s[1]))
			txdatam.append(fm*float(s[1]))
			txdatap.append(fp*float(s[1]))
			tydata.append(float(s[2]))
			tydatam.append(fm*float(s[2]))
			tydatap.append(fp*float(s[2]))
			ttdata.append(float(s[0]))


	m0 = max(ttdata)
	m1 = max(txdata)
	m2 = max(tydata)
	n0 = min(ttdata)
	n1 = min(txdata)
	n2 = min(tydata)

	exdata = []
	ex1data = []
	ex2data = []
	eydata = []
	ey1data = []
	ey2data = []
	etdata = []
	d2 = open("pmpar_e", "r").readlines()
	for d in d2:
		s = split(strip(d))
		if len(s) > 4:
			etdata.append(float(s[0]))
			exdata.append(float(s[1]))
			eydata.append(float(s[3]))
			ex1data.append(float(s[1])-float(s[2]))
			ex2data.append(float(s[1])+float(s[2]))
			ey1data.append(float(s[3])-float(s[4]))
			ey2data.append(float(s[3])+float(s[4]))
	
	m = min(etdata)
	if m < n0 : n0 = m
	m = max(etdata)
	if m > m0 : m0 = m
	m = min(ex1data)
	if m < n1 : n1 = m
	m = max(ex2data)
	if m > m1 : m1 = m
	m = min(ey1data)
	if m < n2 : n2 = m
	m = max(ey2data)
	if m > m2 : m2 = m

	m = m*1.1

	r1 = m1 - n1
	r2 = m2 - n2
	if r1 > r2:
		m2 = m2 + 0.5*(r1-r2)
		n2 = n2 - 0.5*(r1-r2)
	else:
		m1 = m1 + 0.5*(r2-r1)
		n1 = n1 - 0.5*(r2-r1)
	r = m1 - n1
	
	px = biggles.FramedPlot()
	py = biggles.FramedPlot()
	px.yrange = n1-0.05*r, m1+0.05*r
	py.yrange = n2-0.05*r, m2+0.05*r
	px.title = "Motion of %s in RA" % name
	py.title = "Motion of %s in Dec" % name
	px.xlabel = "Date [MJD]"
	px.ylabel = r"$\Delta\alpha$ [mas]"
	py.xlabel = "Date [MJD]"
	py.ylabel = r"$\Delta\delta$ [mas]"
	px.add(biggles.Curve(ttdata, txdata))
	px.add(biggles.YErrorBars(etdata, ex1data, ex2data))
	py.add(biggles.Curve(ttdata, tydata))
	py.add(biggles.YErrorBars(etdata, ey1data, ey2data))
	if showsigpi == 1:
		px.add(biggles.Curve(ttdata, txdatam, linetype='dashed'))
		px.add(biggles.Curve(ttdata, txdatap, linetype='dashed'))
		py.add(biggles.Curve(ttdata, tydatam, linetype='dashed'))
		py.add(biggles.Curve(ttdata, tydatap, linetype='dashed'))
	p = biggles.Table( 2, 1, cellspacing=0.02 )
	p.set(0, 0, px)
	p.set(1, 0, py)
	p.comment = name

	fspace = 0.09
	if find(sys.argv[-1], "S") > 0:
		px.add(biggles.PlotLabel(0.05, 1.0-fspace, mu_a_str, \
			valign = 'top', halign = 'left'))
		px.add(biggles.PlotLabel(0.05, 1.0-2.0*fspace, mu_d_str, \
			valign = 'top', halign = 'left'))
		px.add(biggles.PlotLabel(0.05, 1.0-3.0*fspace, pi_str, \
			valign = 'top', halign = 'left'))

	return p

def multiplotparallax(argv):
	if find(sys.argv[-1], "e") > 0:
		equal = 1
	else:
		equal = 0
	files = argv[1:-1]
	n = round(sqrt(len(files))+0.499999)
	p = biggles.Table(n, n, cellspacing=0.02)
	i = 0
	j = 0
	ax=[]
	bx=[]
	ay=[]
	by=[]
	Plots = []
	for f in files:
		P = makeparallaxplot([argv[0], f, argv[-1]])
		ax.append(P.xrange[0])
		bx.append(P.xrange[1])
		ay.append(P.yrange[0])
		by.append(P.yrange[1])
		Plots.append(P)
	for P in Plots:
		if equal == 1:
			P.xrange = min(ax), max(bx)
			P.yrange = min(ay), max(by)
		p.set(i, j, P)
		i = i+1
		if i >= n:
			i = 0
			j = j+1
	p.comment = "multiplot"
	return p
	
def multiplotparallaxxy(argv):
	files = argv[1:-1]
	n = round(sqrt(len(files))+0.499999)
	p = biggles.Table(n, n, cellspacing=0.02)
	i = 0
	j = 0
	for f in files:
		p.set(i, j, makeparallaxplotxy([argv[0], f, argv[-1]]))
		i = i+1
		if i >= n:
			i = 0
			j = j+1
	p.comment = "multiplot"
	return p

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print "Need more arguments"
	
	if sys.argv[-1][0] != '-':
		sys.argv.append("-o")
	
	if find(sys.argv[-1], "x") > 0:
		doxy = 1
	else:
		doxy = 0
	

	if len(sys.argv) < 4:
		if doxy == 1:
			p = makeparallaxplotxy(sys.argv)
		else:
			p = makeparallaxplot(sys.argv)
	else:
		if doxy == 1:
			p = multiplotparallaxxy(sys.argv)
		else:
			p = multiplotparallax(sys.argv)
	
	if find(sys.argv[-1], "A") > 0:
		print p.xrange
		p.aspect_ratio = 1.618033
		w = p.xrange[1] - p.xrange[0]
		w = w*0.190983
		p.xrange = (p.xrange[0] + w, p.xrange[1] - w)
		autosetcustomticks(p)
		
	if find(sys.argv[-1], "a") > 0:
		p.aspect_ratio = 0.618033
		w = p.yrange[1] - p.yrange[0]
		w = w*0.190983
		p.yrange = (p.yrange[0] + w, p.yrange[1] - w)
		autosetcustomticks(p)

	p.show()
	p.save_as_eps("%s.ps" % p.comment);



