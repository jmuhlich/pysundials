try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math
import sys

NUM_SPECIES = 2
C1_SCALE = 1.0e6
C2_SCALE = 1.0e12

T0 = 0.0
NOUT = 12
TWOHR = 7200.0
HALFDAY = 4.32e4
PI = 3.1415926535898

XMIN = 0.0
XMAX = 20.0
ZMIN = 30.0
ZMAX = 50.0
XMID = 10.0
ZMID = 40.0

MX = 15
MZ = 15
NSMX = NUM_SPECIES*MX
MM = (MX*MZ)

RTOL = 1.0e-5
FLOOR = 100.0
ATOL = (RTOL*FLOOR)
NEQ = (NUM_SPECIES*MM)

NP = 8
NS = 2

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])

class UserData (ctypes.Structure):
	_fields_ = [
		("p", ctypes.POINTER(cvodes.realtype)),
		("P", (ctypes.POINTER(ctypes.POINTER(cvodes.realtype))*MX)*MZ),
		("Jbd", (ctypes.POINTER(ctypes.POINTER(cvodes.realtype))*MX)*MZ),
		("pivot", (ctypes.POINTER(ctypes.c_long)*MX)*MZ),
		("q4", cvodes.realtype),
		("om", cvodes.realtype),
		("dx", cvodes.realtype),
		("dz", cvodes.realtype),
		("hdco", cvodes.realtype),
		("haco", cvodes.realtype),
		("vdco", cvodes.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def SetInitialProfiles(y, dx, dz):
	for jz in range(MZ):
		z = ZMIN + jz*dz
		cz = (0.1*(z - ZMID))**2
		cz = 1.0 - cz + 0.5*(cz**2)
		for jx in range(MX):
			x = XMIN + jx*dx
			cx = (0.1*(x - XMID))**2
			cx = 1.0 - cx + 0.5*cx**2
			y[1-1 + jx*NUM_SPECIES + jz*NSMX] = C1_SCALE*cx*cz
			y[2-1 + jx*NUM_SPECIES + jz*NSMX] = C2_SCALE*cx*cz

def PrintOutput(cvode_mem, t, y):
	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	qu = cvodes.CVodeGetLastOrder(cvode_mem)
	hu = cvodes.CVodeGetLastStep(cvode_mem)
	
	print "%8.3le %2d  %8.3le %5ld"%(t.value,qu,hu,nst)
	
	print "                                Solution      ",
	print "%12.4le %12.4le "%(y[1-1 + (0)*NUM_SPECIES + (0)*NSMX], y[1-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])
	print "                                              ",
	print "%12.4le %12.4le "%(y[2-1 + (0)*NUM_SPECIES + (0)*NSMX], y[2-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])

def PrintOutputS(uS):
	print "                                ----------------------------------------"
	print "                                Sensitivity 1 ",
	print "%12.4le %12.4le "%(uS[0][1-1 + (0)*NUM_SPECIES + (0)*NSMX], uS[0][1-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])
	print "                                              ",
	print "%12.4le %12.4le "%(uS[0][2-1 + (0)*NUM_SPECIES + (0)*NSMX], uS[0][2-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])
	
	print "                                ----------------------------------------"
	print "                                Sensitivity 2 ",
	print "%12.4le %12.4le "%(uS[1][1-1 + (0)*NUM_SPECIES + (0)*NSMX], uS[1][1-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])
	print "                                              ",
	print "%12.4le %12.4le "%(uS[1][2-1 + (0)*NUM_SPECIES + (0)*NSMX], uS[1][2-1 + (MX-1)*NUM_SPECIES + (MZ-1)*NSMX])
	

def PrintFinalStats(cvode_mem, sensi):
	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	nfe = cvodes.CVodeGetNumRhsEvals(cvode_mem)
	nsetups = cvodes.CVodeGetNumLinSolvSetups(cvode_mem)
	netf = cvodes.CVodeGetNumErrTestFails(cvode_mem)
	nni = cvodes.CVodeGetNumNonlinSolvIters(cvode_mem)
	ncfn = cvodes.CVodeGetNumNonlinSolvConvFails(cvode_mem)
	
	if sensi:
		nfSe = cvodes.CVodeGetNumSensRhsEvals(cvode_mem)
		nfeS = cvodes.CVodeGetNumRhsEvalsSens(cvode_mem)
		nsetupsS = cvodes.CVodeGetNumSensLinSolvSetups(cvode_mem)
		netfS = cvodes.CVodeGetNumSensErrTestFails(cvode_mem)
		nniS = cvodes.CVodeGetNumSensNonlinSolvIters(cvode_mem)
		ncfnS = cvodes.CVodeGetNumSensNonlinSolvConvFails(cvode_mem)
	
	nli = cvodes.CVSpilsGetNumLinIters(cvode_mem)
	ncfl = cvodes.CVSpilsGetNumConvFails(cvode_mem)
	npe = cvodes.CVSpilsGetNumPrecEvals(cvode_mem)
	nps = cvodes.CVSpilsGetNumPrecSolves(cvode_mem)
	
	print "\nFinal Statistics\n"
	print "nst     = %5ld\n"%(nst)
	print "nfe     = %5ld"%(nfe)
	print "netf    = %5ld    nsetups  = %5ld"%(netf, nsetups)
	print "nni     = %5ld    ncfn     = %5ld"%(nni, ncfn)
	
	if sensi:
		print
		print "nfSe    = %5ld    nfeS     = %5ld"%(nfSe, nfeS)
		print "netfs   = %5ld    nsetupsS = %5ld"%(netfS, nsetupsS)
		print "nniS    = %5ld    ncfnS    = %5ld"%(nniS, ncfnS)
		
	print
	print "nli     = %5ld    ncfl     = %5ld"%(nli, ncfl)
	print "npe     = %5ld    nps      = %5ld"%(npe, nps)

def f(t, y, ydot, f_data):
	data = ctypes.cast(f_data, PUserData).contents

	#Set diurnal rate coefficients

	s = math.sin(data.om*t)
	if s > 0:
		q3 = math.exp(-data.p[3]/s)
		data.q4 = math.exp(-data.p[4]/s)
	else:
		q3 = 0
		data.q4 = 0

	#Make local copies of problem variables, for efficiency.

	q4coef = data.q4
	delz = data.dz
	verdco = data.vdco
	hordco = data.hdco
	horaco = data.haco

	#Loop over all grid points.

	for jz in range(MZ):

		#Set vertical diffusion coefficients at jz +- 1/2

		zdn = ZMIN + (jz - 0.5)*delz
		zup = zdn + delz
		czdn = verdco*math.exp(0.2*zdn)
		czup = verdco*math.exp(0.2*zup)
		if jz == 0:
			idn = 1
		else:
			idn = -1
		if jz == MZ-1:
			iup = -1
		else:
			iup = 1

		for jx in range(MX):

			#Extract c1 and c2, and set kinetic rate terms.

			c1 = y[1-1 + (jx)*NUM_SPECIES + (jz)*NSMX] 
			c2 = y[2-1 + (jx)*NUM_SPECIES + (jz)*NSMX]
			qq1 = data.p[0]*c1*data.p[2]
			qq2 = data.p[1]*c1*c2
			qq3 = q3*data.p[2]
			qq4 = q4coef*c2
			rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4
			rkin2 = qq1 - qq2 - qq4

			#Set vertical diffusion terms.

			c1dn = y[1-1 + (jx)*NUM_SPECIES + (jz+idn)*NSMX]
			c2dn = y[2-1 + (jx)*NUM_SPECIES + (jz+idn)*NSMX]
			c1up = y[1-1 + (jx)*NUM_SPECIES + (jz+iup)*NSMX]
			c2up = y[2-1 + (jx)*NUM_SPECIES + (jz+iup)*NSMX]
			vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn)
			vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn)

			#Set horizontal diffusion and advection terms.

			if jx == 0:
				ileft = 1
			else:
				ileft = -1
			if jx == MX-1:
				iright = -1
			else:
				iright = 1
			c1lt = y[1-1 + (jx+ileft)*NUM_SPECIES + (jz)*NSMX]
			c2lt = y[2-1 + (jx+ileft)*NUM_SPECIES + (jz)*NSMX]
			c1rt = y[1-1 + (jx+iright)*NUM_SPECIES + (jz)*NSMX]
			c2rt = y[2-1 + (jx+iright)*NUM_SPECIES + (jz)*NSMX]
			hord1 = hordco*(c1rt - 2.0*c1 + c1lt)
			hord2 = hordco*(c2rt - 2.0*c2 + c2lt)
			horad1 = horaco*(c1rt - c1lt)
			horad2 = horaco*(c2rt - c2lt)

			#Load all terms into udot.

			ydot[1-1 + (jx)*NUM_SPECIES + (jz)*NSMX] = vertd1 + hord1 + horad1 + rkin1
			ydot[2-1 + (jx)*NUM_SPECIES + (jz)*NSMX] = vertd2 + hord2 + horad2 + rkin2

	return 0

def Precond(tn, y, fy, jok, jcurPtr, gamma, P_data, vtemp1, vtemp2, vtemp3):
	data = ctypes.cast(P_data, PUserData).contents
	
	if jok:
		for jz in range(MZ):
			for jx in range(MZ):
				cvodes.dencopy(data.Jbd[jx][jz], data.P[jx][jz], NUM_SPECIES, NUM_SPECIES)
		jcurPtr.contents.value = 0
	else:
		q4coef = data.q4
		delz = data.dz
		verdco = data.vdco
		hordco = data.hdco
		
		for jz in range(MZ):
			zdn = ZMIN + (jz - 0.5)*delz
			zup = zdn + delz
			czdn = verdco*math.exp(0.2*zdn)
			czup = verdco*math.exp(0.2*zup)
			diag = -(czdn + czup + 2.0*hordco)
			for jx in range(MX):
				c1 = y[1-1 + (jx)*NUM_SPECIES + (jz)*NSMX]
				c2 = y[2-1 + (jx)*NUM_SPECIES + (jz)*NSMX]
				j = data.Jbd[jx][jz]
				a = data.P[jx][jz]
				j[0][0] = (-data.p[0]*data.p[2] - data.p[1]*c2) + diag
				j[1][0] = -data.p[1]*c1 + q4coef
				j[0][1] = data.p[0]*data.p[2] - data.p[1]*c2
				j[1][1] = (-data.p[1]*c1 - q4coef) + diag
				cvodes.dencopy(j, a, NUM_SPECIES, NUM_SPECIES)
		jcurPtr.contents.value = 1
	
	for jz in range(MZ):
		for jx in range(MX):
			cvodes.denscale(-gamma, data.P[jx][jz], NUM_SPECIES, NUM_SPECIES)
	
	for jx in range(MX):
		for jz in range(MZ):
			cvodes.denaddI(data.P[jx][jz], NUM_SPECIES)
			ier = cvodes.denGETRF(data.P[jx][jz], NUM_SPECIES, NUM_SPECIES, data.pivot[jx][jz])
			if ier != 0:
				return 1
	
	return 0

def PSolve(tn, y, fy, r, z, gamma, delta, lr, P_data, vtemp):
	data = ctypes.cast(P_data, PUserData).contents
	
	z[:] = r

	for jx in range(MX):
		for jz in range(MZ):
			cvodes.denGETRS(data.P[jx][jz], NUM_SPECIES, data.pivot[jx][jz], z.ptrto(jx*NUM_SPECIES + jz*NSMX))

	return 0

def WrongArgs(name):
	print "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]"%(name)
	print "         sensi_meth = sim, stg, or stg1"
	print "         err_con    = t or f"

	sys.exit(0)

sensi = False
sensi_meth = -1
err_con = False
t = cvodes.realtype(0)

if len(sys.argv) < 2:
	WrongArgs(sys.argv[0])

if sys.argv[1] == "-nosensi":
	sensi = False
elif sys.argv[1] == "-sensi":
	sensi = True
else:
	WrongArgs(sys.argv[0])

if (sensi):
	if len(sys.argv) != 4:
		WrongArgs(sys.argv[0])

	if sys.argv[2] == "sim":
		sensi_meth = cvodes.CV_SIMULTANEOUS
	elif sys.argv[2] == "stg":
		sensi_meth = cvodes.CV_STAGGERED
	elif sys.argv[2] == "stg1":
		sensi_meth = cvodes.CV_STAGGERED1
	else:
		WrongArgs(sys.argv[0])

	if sys.argv[3] == "t":
		err_con = True
	elif sys.argv[3] == "f":
		err_con = False
	else:
		WrongArgs(sys.argv[0])

#Allocate and initialise user data
t = cvodes.realtype(0)
data = UserData()
for jx in range(MX):
	for jz in range(MZ):
		data.P[jx][jz] = cvodes.denalloc(NUM_SPECIES, NUM_SPECIES)
		data.Jbd[jx][jz] = cvodes.denalloc(NUM_SPECIES, NUM_SPECIES)
		data.pivot[jx][jz] = cvodes.denallocpiv(NUM_SPECIES)
params = (cvodes.realtype*NP)()

data.om = PI/HALFDAY
data.dx = (XMAX-XMIN)/(MX-1)
data.dz = (ZMAX-ZMIN)/(MZ-1)
data.hdco = 4.0e-6/(data.dx**2)
data.haco = 0.001/(2.0*data.dx)
data.vdco = (1.0/(data.dz**2))*1.0e-8
data.p = params

data.p[0] = 1.63e-16
data.p[1] = 4.66e-16
data.p[2] = 3.7e16
data.p[3] = 22.62
data.p[4] = 7.601
data.p[5] = 4.0e-6
data.p[6] = 0.001
data.p[7] = 1.0e-8

y = cvodes.NVector([0]*(NEQ))
SetInitialProfiles(y, data.dx, data.dz)
  
cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)

cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVodeSetMaxNumSteps(cvode_mem, 2000)
cvodes.CVodeMalloc(cvode_mem, f, T0, y, cvodes.CV_SS, RTOL, ATOL)
cvodes.CVSpgmr(cvode_mem, cvodes.PREC_LEFT, 0)
cvodes.CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve, ctypes.pointer(data))

print "\n2-species diurnal advection-diffusion problem"

if sensi:
	plist = (ctypes.c_int*NS)()
	for i in range(NS):
		plist[i] = i

	pbar = (cvodes.realtype*NS)(data.p[0], data.p[1])

	uS = nvecserial.NVectorArray([[0]*NEQ]*NS)

	cvodes.CVodeSensMalloc(cvode_mem, NS, sensi_meth, uS)
	cvodes.CVodeSetSensErrCon(cvode_mem, err_con)
	cvodes.CVodeSetSensDQMethod(cvode_mem, cvodes.CV_CENTERED, 0)
	cvodes.CVodeSetSensParams(cvode_mem, data.p, pbar, plist)

	print "Sensitivity: YES",
	if sensi_meth == cvodes.CV_SIMULTANEOUS:
		print "( SIMULTANEOUS +",
	else:
		if sensi_meth == cvodes.CV_STAGGERED:
			print "( STAGGERED +", 
		else:
			print "( STAGGERED1 +",
	if err_con:
		print "FULL ERROR CONTROL )"
	else:
		print "PARTIAL ERROR CONTROL )"

else:
	print "Sensitivity: NO "
	print "\n"
	print "========================================================================"
	print "     T     Q       H      NST                    Bottom left  Top right "
	print "========================================================================"

iout = 1
tout = TWOHR
while iout <= NOUT:
	cvodes.CVode(cvode_mem, tout, y, ctypes.byref(t), cvodes.CV_NORMAL)
	PrintOutput(cvode_mem, t, y)

	if sensi:
		cvodes.CVodeGetSens(cvode_mem, t, uS)
		PrintOutputS(uS)
	print "------------------------------------------------------------------------"

	iout += 1
	tout += TWOHR

PrintFinalStats(cvode_mem, sensi)

for jx in range(MX):
	for jz in range(MZ):
		cvodes.denfree(data.P[jx][jz])
		cvodes.denfree(data.Jbd[jx][jz])
		cvodes.denfreepiv(data.pivot[jx][jz])
