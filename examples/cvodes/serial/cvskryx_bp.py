try:
	from pysundials import cvodes
except ImportError:
	import cvodes
import ctypes
import math

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

T0 = 0
NOUT = 12
TWOHR = 7200.0
HALFDAY = 4.32e4

XMIN = 0
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

class UserData(ctypes.Structure):
	_fields_ = [
		('q4', cvodes.realtype),
		('om', cvodes.realtype),
		('dx', cvodes.realtype),
		('dy', cvodes.realtype),
		('hdco', cvodes.realtype),
		('haco', cvodes.realtype),
		('vdco', cvodes.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def InitUserData(data):
	data.om = math.pi/HALFDAY
	data.dx = (XMAX-XMIN)/(MX-1)
	data.dy = (YMAX-YMIN)/(MY-1)
	data.hdco = KH/(data.dx)**2
	data.haco = VEL/(1.0*data.dx)
	data.vdco = (1.0/(data.dy)**2)*KV0

def SetInitialProfiles(u, dx, dy):
	for jy in range(MY):
		y = YMIN + jy*dy
		cy = (0.1*(y - YMID))**2
		cy = 1.0 - cy + 0.5*(cy)**2
		for jx in range(MX):
			x = XMIN + jx*dx
			cx = (0.1*(x - XMID))**2
			cx = 1.0 - cx + 0.5*(cx)**2
			u[1-1 + (jx)*NUM_SPECIES + (jy)*NSMX] = C1_SCALE*cx*cy
			u[2-1 + (jx)*NUM_SPECIES + (jy)*NSMX] = C2_SCALE*cx*cy

def PrintIntro(mu, ml):
	print "2-species diurnal advection-diffusion problem, %d by %d mesh"%(MX, MY)
	print "SPGMR solver; band preconditioner; mu = %d, ml = %d\n"%(mu, ml)

def PrintOutput(cvode_mem, u, t):
	mxh = MX/2 - 1 
	myh = MY/2 - 1
	mx1 = MX - 1
	my1 = MY - 1

	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	qu = cvodes.CVodeGetLastOrder(cvode_mem)
	hu = cvodes.CVodeGetLastStep(cvode_mem)

	print "t = %.2le   no. steps = %ld  order = %d  stepsize = %.2le"%(t.value, nst, qu, hu)
	print "c1 (bot.left/middle/top rt.) = %12.3le  %12.3le  %12.3le"%(u[1-1 + (0)*NUM_SPECIES + (0)*NSMX], u[1-1 + (mxh)*NUM_SPECIES + (myh)*NSMX], u[1-1 + (mx1)*NUM_SPECIES + (my1)*NSMX])
	print "c2 (bot.left/middle/top rt.) = %12.3le  %12.3le  %12.3le\n"%(u[2-1 + (0)*NUM_SPECIES + (0)*NSMX], u[2-1 + (mxh)*NUM_SPECIES + (myh)*NSMX], u[2-1 + (mx1)*NUM_SPECIES + (my1)*NSMX])

def PrintFinalStats(cvode_mem, bpdata):
	(lenrw, leniw) = cvodes.CVodeGetWorkSpace(cvode_mem)
	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	nfe = cvodes.CVodeGetNumRhsEvals(cvode_mem)
	nsetups = cvodes.CVodeGetNumLinSolvSetups(cvode_mem)
	netf = cvodes.CVodeGetNumErrTestFails(cvode_mem)
	nni = cvodes.CVodeGetNumNonlinSolvIters(cvode_mem)
	ncfn = cvodes.CVodeGetNumNonlinSolvConvFails(cvode_mem)

	(lenrwLS, leniwLS) = cvodes.CVSpilsGetWorkSpace(cvode_mem)
	nli = cvodes.CVSpilsGetNumLinIters(cvode_mem)
	npe = cvodes.CVSpilsGetNumPrecEvals(cvode_mem)
	nps = cvodes.CVSpilsGetNumPrecSolves(cvode_mem)
	ncfl = cvodes.CVSpilsGetNumConvFails(cvode_mem)
	nfeLS = cvodes.CVSpilsGetNumRhsEvals(cvode_mem)

	(lenrwBP, leniwBP) = cvodes.CVBandPrecGetWorkSpace(bpdata)
	nfeBP = cvodes.CVBandPrecGetNumRhsEvals(bpdata)

	print "\nFinal Statistics.. \n"
	print "lenrw   = %5ld     leniw   = %5ld"%(lenrw.value, leniw.value)
	print "lenrwls = %5ld     leniwls = %5ld"%(lenrwLS, leniwLS)
	print "lenrwbp = %5ld     leniwbp = %5ld"%(lenrwBP, leniwBP)
	print "nst     = %5ld"%(nst)
	print "nfe     = %5ld     nfetot  = %5ld"%(nfe, nfe+nfeLS+nfeBP)
	print "nfeLS   = %5ld     nfeBP   = %5ld"%(nfeLS, nfeBP)
	print "nni     = %5ld     nli     = %5ld"%(nni, nli)
	print "nsetups = %5ld     netf    = %5ld"%(nsetups, netf)
	print "npe     = %5ld     nps     = %5ld"%(npe, nps)
	print "ncfn    = %5ld     ncfl    = %5ld\n"%(ncfn, ncfl)

def f(t, u, udot, f_data):
	data = ctypes.cast(f_data, PUserData).contents

	s = math.sin(data.om*t)
	if s > 0:
		q3 = math.exp(-A3/s)
		data.q4 = math.exp(-A4/s)
	else:
		q3 = 0
		data.q4 = 0

	q4coef = data.q4
	dely = data.dy
	verdco = data.vdco
	hordco = data.hdco
	horaco = data.haco

	for jy in range(MY):
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
			c1 = u[1-1 + (jx)*NUM_SPECIES + (jy)*NSMX]
			c2 = u[2-1 + (jx)*NUM_SPECIES + (jy)*NSMX]
			qq1 = Q1*c1*C3
			qq2 = Q2*c1*c2
			qq3 = q3*C3
			qq4 = q4coef*c2
			rkin1 = -qq1 - qq2 + 1.0*qq3 + qq4
			rkin2 = qq1 - qq2 - qq4

			c1dn = u[1-1 + (jx)*NUM_SPECIES + (jy+idn)*NSMX]
			c2dn = u[2-1 + (jx)*NUM_SPECIES + (jy+idn)*NSMX]
			c1up = u[1-1 + (jx)*NUM_SPECIES + (jy+iup)*NSMX]
			c2up = u[2-1 + (jx)*NUM_SPECIES + (jy+iup)*NSMX]
			vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn)
			vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn)

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
			hord1 = hordco*(c1rt - 1.0*c1 + c1lt)
			hord2 = hordco*(c2rt - 1.0*c2 + c2lt)
			horad1 = horaco*(c1rt - c1lt)
			horad2 = horaco*(c2rt - c2lt)

			udot[ 1-1 + ( jx)*NUM_SPECIES + ( jy)*NSMX] = vertd1 + hord1 + horad1 + rkin1
			udot[ 2-1 + ( jx)*NUM_SPECIES + ( jy)*NSMX] = vertd2 + hord2 + horad2 + rkin2
	return 0

t = cvodes.realtype(0)
u = cvodes.NVector([0]*NEQ)
data = UserData()
InitUserData(data)
SetInitialProfiles(u, data.dx, data.dy)
abstol = cvodes.realtype(ATOL)
reltol = RTOL

cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVodeMalloc(cvode_mem, f, T0, u, cvodes.CV_SS, reltol, abstol)
ml = mu = 2
bpdata = cvodes.CVBandPrecAlloc(cvode_mem, NEQ, mu, ml)
cvodes.CVBPSpgmr(cvode_mem, cvodes.PREC_LEFT, 0, bpdata)
PrintIntro(mu, ml)

for jpre in range(cvodes.PREC_LEFT, cvodes.PREC_RIGHT+1):
	if jpre == cvodes.PREC_RIGHT:
		SetInitialProfiles(u, data.dx, data.dy)
		cvodes.CVodeReInit(cvode_mem, f, T0, u, cvodes.CV_SS, reltol, abstol)
		cvodes.CVSpilsSetPrecType(cvode_mem, cvodes.PREC_RIGHT)
		print "\n\n-------------------------------------------------------------------"
	
		print "\n\nPreconditioner type is:  jpre = PREC_RIGHT\n"
	else:
		print "\n\nPreconditioner type is:  jpre = PREC_LEFT\n"
	
	iout = 1
	tout = TWOHR
	while iout <= NOUT:
		flag = cvodes.CVode(cvode_mem, tout, u, ctypes.byref(t), cvodes.CV_NORMAL)
		PrintOutput(cvode_mem, u, t)
		if flag != cvodes.CV_SUCCESS:
			break
		iout += 1
		tout += TWOHR
	PrintFinalStats(cvode_mem, bpdata)
