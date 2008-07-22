#!/usr/bin/python

try:
	from pysundials import cvode
except ImportError:
	import cvode
import math
import ctypes

NUM_SPECIES = 2
KH = 4.0e-6
VEL = 0.001
KV0 = 1.0e-8
Q1 = 1.63e-16
Q2 = 4.66e-16
C3 = 3.7e16
A3 = 22.62
A4 = 7.601
C1_SCALE = 1.0e6
C2_SCALE = 1.0e12

T0 = 0.0
NOUT = 12
TWOHR = 7200.0
HALFDAY = 4.32e4
PI = 3.1415926535898

XMIN = 0.0
XMAX = 20.0
YMIN = 30.0
YMAX = 50.0
XMID = 10.0
YMID = 40.0

MX = 10
MY = 10
NSMX = 20
MM = (MX*MY)

RTOL = 1.0e-5
FLOOR = 100.0
ATOL = (RTOL*FLOOR)
NEQ = (NUM_SPECIES*MM)

USE_SPGMR = 0
USE_SPBCG = 1
USE_SPTFQMR = 2

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])

class UserData (ctypes.Structure):
	_fields_ = [
		("P", (ctypes.POINTER(ctypes.POINTER(cvode.realtype))*MX)*MY),
		("Jbd", (ctypes.POINTER(ctypes.POINTER(cvode.realtype))*MX)*MY),
		("pivot", (ctypes.POINTER(ctypes.c_long)*MX)*MY),
		("q4", cvode.realtype),
		("om", cvode.realtype),
		("dx", cvode.realtype),
		("dy", cvode.realtype),
		("hdco", cvode.realtype),
		("haco", cvode.realtype),
		("vdco", cvode.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def SetInitialProfiles(u, dx, dy):
	for jy in range(MY):
		y = YMIN + jy*dy
		cy = (0.1*(y - YMID))**2
		cy = 1.0 - cy + 0.5*(cy**2)
		for jx in range(MX):
			x = XMIN + jx*dx
			cx = (0.1*(x - XMID))**2
			cx = 1.0 - cx + 0.5*cx**2
			u[1-1 + jx*NUM_SPECIES + jy*NSMX] = C1_SCALE*cx*cy
			u[2-1 + jx*NUM_SPECIES + jy*NSMX] = C2_SCALE*cx*cy

def PrintOutput(cvode_mem, u, t):
	mxh = MX/2 - 1
	myh = MY/2 - 1
	mx1 = MX - 1
	my1 = MY - 1

	nst = cvode.CVodeGetNumSteps(cvode_mem)
	qu = cvode.CVodeGetLastOrder(cvode_mem)
	hu = cvode.CVodeGetLastStep(cvode_mem)

	print "t = %8.2e   no. steps = %i   order = %i   stepsize = %.2g"%(t.value, nst, qu, hu)
	print "c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e"%(u[1-1 + (0)*NUM_SPECIES + (0)*NSMX], u[1-1 + (mxh)*NUM_SPECIES + (myh)*NSMX], u[1-1 + (mx1)*NUM_SPECIES + (my1)*NSMX])
	print "c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n"%(u[2-1 + (0)*NUM_SPECIES + (0)*NSMX], u[2-1 + (mxh)*NUM_SPECIES + (myh)*NSMX], u[2-1 + (mx1)*NUM_SPECIES + (my1)*NSMX])

def PrintFinalStats(cvode_mem, linsolver):
	(lenrw, leniw) = cvode.CVodeGetWorkSpace(cvode_mem)
	nst = cvode.CVodeGetNumSteps(cvode_mem)
	nfe = cvode.CVodeGetNumRhsEvals(cvode_mem)
	nsetups = cvode.CVodeGetNumLinSolvSetups(cvode_mem)
	netf = cvode.CVodeGetNumErrTestFails(cvode_mem)
	nni = cvode.CVodeGetNumNonlinSolvIters(cvode_mem)
	ncfn = cvode.CVodeGetNumNonlinSolvConvFails(cvode_mem)
	
	(lenrwLS, leniwLS) = cvode.CVSpilsGetWorkSpace(cvode_mem)
	nli = cvode.CVSpilsGetNumLinIters(cvode_mem)
	npe = cvode.CVSpilsGetNumPrecEvals(cvode_mem)
	nps = cvode.CVSpilsGetNumPrecSolves(cvode_mem)
	ncfl = cvode.CVSpilsGetNumConvFails(cvode_mem)
	nfeLS = cvode.CVSpilsGetNumRhsEvals(cvode_mem)
	
	print "\nFinal Statistics.. \n"
	print "lenrw   = %5ld     leniw   = %5ld"%(lenrw.value, leniw.value)
	print "lenrwLS = %5ld     leniwLS = %5ld"%(lenrwLS, leniwLS)
	print "nst     = %5ld"%(nst)
	print "nfe     = %5ld     nfeLS   = %5ld"%(nfe, nfeLS)
	print "nni     = %5ld     nli     = %5ld"%(nni, nli)
	print "nsetups = %5ld     netf    = %5ld"%(nsetups, netf)
	print "npe     = %5ld     nps     = %5ld"%(npe, nps)
	print "ncfn    = %5ld     ncfl    = %5ld\n"%(ncfn, ncfl)
  
  	if linsolver < 2:
	    print "======================================================================\n"

def f(t, u, udot, f_data):
	data = ctypes.cast(f_data, PUserData).contents

	#Set diurnal rate coefficients

	s = math.sin(data.om*t)
	if s > 0:
		q3 = math.exp(-A3/s)
		data.q4 = math.exp(-A4/s)
	else:
		q3 = 0
		data.q4 = 0

	#Make local copies of problem variables, for efficiency.

	q4coef = data.q4
	dely = data.dy
	verdco = data.vdco
	hordco = data.hdco
	horaco = data.haco

	#Loop over all grid points.

	for jy in range(MY):

		#Set vertical diffusion coefficients at jy +- 1/2

		ydn = YMIN + (jy - 0.5)*dely
		yup = ydn + dely
		cydn = verdco*math.exp(0.2*ydn)
		cyup = verdco*math.exp(0.2*yup)
		if jy == 0:
			idn = 1
		else:
			idn = -1
		if jy == MY-1:
			iup = -1
		else:
			iup = 1

		for jx in range(MX):

			#Extract c1 and c2, and set kinetic rate terms.

			c1 = u[1-1 + (jx)*NUM_SPECIES + (jy)*NSMX] 
			c2 = u[2-1 + (jx)*NUM_SPECIES + (jy)*NSMX]
			qq1 = Q1*c1*C3
			qq2 = Q2*c1*c2
			qq3 = q3*C3
			qq4 = q4coef*c2
			rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4
			rkin2 = qq1 - qq2 - qq4

			#Set vertical diffusion terms.

			c1dn = u[1-1 + (jx)*NUM_SPECIES + (jy+idn)*NSMX]
			c2dn = u[2-1 + (jx)*NUM_SPECIES + (jy+idn)*NSMX]
			c1up = u[1-1 + (jx)*NUM_SPECIES + (jy+iup)*NSMX]
			c2up = u[2-1 + (jx)*NUM_SPECIES + (jy+iup)*NSMX]
			vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn)
			vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn)

			#Set horizontal diffusion and advection terms.

			if jx == 0:
				ileft = 1
			else:
				ileft = -1
			if jx == MX-1:
				iright = -1
			else:
				iright = 1
			c1lt = u[1-1 + (jx+ileft)*NUM_SPECIES + (jy)*NSMX]
			c2lt = u[2-1 + (jx+ileft)*NUM_SPECIES + (jy)*NSMX]
			c1rt = u[1-1 + (jx+iright)*NUM_SPECIES + (jy)*NSMX]
			c2rt = u[2-1 + (jx+iright)*NUM_SPECIES + (jy)*NSMX]
			hord1 = hordco*(c1rt - 2.0*c1 + c1lt)
			hord2 = hordco*(c2rt - 2.0*c2 + c2lt)
			horad1 = horaco*(c1rt - c1lt)
			horad2 = horaco*(c2rt - c2lt)

			#Load all terms into udot.

			udot[1-1 + (jx)*NUM_SPECIES + (jy)*NSMX] = vertd1 + hord1 + horad1 + rkin1
			udot[2-1 + (jx)*NUM_SPECIES + (jy)*NSMX] = vertd2 + hord2 + horad2 + rkin2

	return 0

def Precond(tn, u, fu, jok, jcurPtr, gamma, P_data, vtemp1, vtemp2, vtemp3):
	data = ctypes.cast(P_data, PUserData).contents
	
	if jok:
		for jy in range(MY):
			for jx in range(MY):
				cvode.dencopy(data.Jbd[jx][jy], data.P[jx][jy], NUM_SPECIES, NUM_SPECIES)
		jcurPtr.contents.value = 0
	else:
		q4coef = data.q4
		dely = data.dy
		verdco = data.vdco
		hordco = data.hdco
		
		for jy in range(MY):
			ydn = YMIN + (jy - 0.5)*dely
			yup = ydn + dely
			cydn = verdco*math.exp(0.2*ydn)
			cyup = verdco*math.exp(0.2*yup)
			diag = -(cydn + cyup + 2.0*hordco)
			for jx in range(MX):
				c1 = u[1-1 + (jx)*NUM_SPECIES + (jy)*NSMX]
				c2 = u[2-1 + (jx)*NUM_SPECIES + (jy)*NSMX]
				j = data.Jbd[jx][jy]
				a = data.P[jx][jy]
				j[0][0] = (-Q1*C3 - Q2*c2) + diag
				j[1][0] = -Q2*c1 + q4coef
				j[0][1] = Q1*C3 - Q2*c2
				j[1][1] = (-Q2*c1 - q4coef) + diag
				cvode.dencopy(j, a, NUM_SPECIES, NUM_SPECIES)
		jcurPtr.contents.value = 1
	
	for jy in range(MY):
		for jx in range(MX):
			cvode.denscale(-gamma, data.P[jx][jy], NUM_SPECIES, NUM_SPECIES)
	
	for jx in range(MX):
		for jy in range(MY):
			cvode.denaddI(data.P[jx][jy], NUM_SPECIES)
			ier = cvode.denGETRF(data.P[jx][jy], NUM_SPECIES, NUM_SPECIES, data.pivot[jx][jy])
			if ier != 0:
				return 1
	
	return 0

def PSolve(tn, u, fu, r, z, gamma, delta, lr, P_data, vtemp):
	data = ctypes.cast(P_data, PUserData).contents
	
	z[:] = r

	for jx in range(MX):
		for jy in range(MY):
			cvode.denGETRS(data.P[jx][jy], NUM_SPECIES, data.pivot[jx][jy], z.ptrto(jx*NUM_SPECIES + jy*NSMX))
	return 0

u = cvode.NVector([0.0]*(NEQ))

#Allocate and initialise user data
t = cvode.realtype(0)
data = UserData()
for jx in range(MX):
	for jy in range(MY):
		data.P[jx][jy] = cvode.denalloc(NUM_SPECIES, NUM_SPECIES)
		data.Jbd[jx][jy] = cvode.denalloc(NUM_SPECIES, NUM_SPECIES)
		data.pivot[jx][jy] = cvode.denallocpiv(NUM_SPECIES)

data.om = PI/HALFDAY
data.dx = (XMAX-XMIN)/(MX-1)
data.dy = (YMAX-YMIN)/(MY-1)
data.hdco = KH/(data.dx**2)
data.haco = VEL/(2.0*data.dx)
data.vdco = (1.0/(data.dy**2))*KV0
pdata = ctypes.pointer(data)

SetInitialProfiles(u, data.dx, data.dy)
abstol= cvode.realtype(1.0e-5*100)
reltol= cvode.realtype(1.0e-5)

cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
cvode.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvode.CVodeMalloc(cvode_mem, f, T0, u, cvode.CV_SS, reltol, abstol)

for linsolver in range(3):
	if linsolver != 0:
		data.om = PI/HALFDAY
		data.dx = (XMAX-XMIN)/(MX-1)
		data.dy = (YMAX-YMIN)/(MY-1)
		data.hdco = KH/(data.dx**2)
		data.haco = VEL/(2.0*data.dx)
		data.vdco = (1.0/(data.dy**2))*KV0
		pdata = ctypes.pointer(data)

		SetInitialProfiles(u, data.dx, data.dy)
		cvode.CVodeReInit(cvode_mem, f, T0, u, cvode.CV_SS, reltol, abstol)

	if linsolver == USE_SPGMR:
		print " ------- "
		print "| SPGMR |"
		print " -------"
		cvode.CVSpgmr(cvode_mem, cvode.PREC_LEFT, 0)
		cvode.CVSpilsSetGSType(cvode_mem, cvode.MODIFIED_GS)
		cvode.CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve, ctypes.pointer(data))
	elif  linsolver == USE_SPBCG:
		print " ------- "
		print "| SPBCG |"
		print " -------"
		cvode.CVSpbcg(cvode_mem, cvode.PREC_LEFT, 0)
		cvode.CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve, ctypes.pointer(data))
	elif linsolver == USE_SPTFQMR:
		print " --------- "
		print "| SPTFQMR |"
		print " ---------"
		cvode.CVSptfqmr(cvode_mem, cvode.PREC_LEFT, 0)
		cvode.CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve, ctypes.pointer(data))

	print " \n2-species diurnal advection-diffusion problem\n"
	iout = 1
	tout = TWOHR
	while iout <= NOUT:
		cvode.CVode(cvode_mem, tout, u, ctypes.byref(t), cvode.CV_NORMAL)
		PrintOutput(cvode_mem, u, t)
		iout += 1
		tout += TWOHR
	PrintFinalStats(cvode_mem, linsolver)

for jx in range(MX):
	for jy in range(MY):
		cvode.denfree(data.P[jx][jy])
		cvode.denfree(data.Jbd[jx][jy])
		cvode.denfreepiv(data.pivot[jx][jy])
