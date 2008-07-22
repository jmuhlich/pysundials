try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math
import sys

NEQ = 3
Y1 = 1.0
Y2 = 0.0
Y3 = 0.0
RTOL = 1e-4
ATOL1 = 1e-8
ATOL2 = 1e-14
ATOL3 = 1e-6
T0 = 0.0
T1 = 0.4
TMULT = 10.0
NOUT = 12

NP = 3
NS = 3

class UserData(ctypes.Structure):
	_fields_ = [
		('p', cvodes.realtype*3)
	]
PUserData = ctypes.POINTER(UserData)

def f(t, y, ydot, f_data):
	data = ctypes.cast(f_data, PUserData).contents

	ydot[0] = -data.p[0]*y[0] + data.p[1]*y[1]*y[2]
	ydot[2] = data.p[2]*y[1]*y[1]
	ydot[1] = -ydot[0] - ydot[2]

	return 0

def Jac(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):
	data = ctypes.cast(jac_data, PUserData).contents

	J[0][0] = -data.p[0]
	J[0][1] = data.p[1]*y[2]
	J[0][2] = data.p[1]*y[1]
	J[1][0] =	data.p[0]
	J[1][1] = -data.p[1]*y[2]-2*data.p[2]*y[1]
	J[1][2] = -data.p[1]*y[1]
	J[2][1] = 2*data.p[2]*y[1]
	return 0

def fS(Ns, t, y, ydot, iS, yS, ySdot, fS_data, tmp1, tmp2):
	data = ctypes.cast(fS_data, PUserData).contents

	sd1 = -data.p[0]*yS[0] + data.p[1]*y[2]*yS[1] + data.p[1]*y[1]*yS[2]
	sd3 = 2*data.p[2]*y[1]*yS[1]
	sd2 = -sd1-sd3

	if iS == 0:
		sd1 += -y[0]
		sd2 +=	y[0]
	elif iS == 1:
		sd1 +=	y[1]*y[2]
		sd2 += -y[1]*y[2]
	elif iS == 2:
		sd2 += -y[1]*y[1]
		sd3 +=	y[1]*y[1]

	ySdot[0] = sd1
	ySdot[1] = sd2
	ySdot[2] = sd3

	return 0

def ewt(y, w, e_data):
	atol = [ATOL1, ATOL2, ATOL3]

	for i in range(3):
		ww = RTOL * abs(y[i]) + atol[i];
		if (ww <= 0.0):
			return -1
		w[i] = 1.0/ww

	return 0

def WrongArgs(name):
	print "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]"%(name)
	print "         sensi_meth = sim, stg, or stg1"
	print "         err_con    = t or f"

	sys.exit(0)

def PrintOutput(cvode_mem, t, u):
	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	qu = cvodes.CVodeGetLastOrder(cvode_mem)
	hu = cvodes.CVodeGetLastStep(cvode_mem)

	print "%8.3le %2d  %8.3le %5ld"%(t.value, qu, hu, nst)
	print "                  Solution      ",
	print "%12.4le %12.4le %12.4le "%(u[0], u[1], u[2])

def PrintOutputS(uS):
	print "                  Sensitivity 1  %12.4le %12.4le %12.4le "%(uS[0][0], uS[0][1], uS[0][2])
	print "                  Sensitivity 2  %12.4le %12.4le %12.4le "%(uS[1][0], uS[1][1], uS[1][2])
	print "                  Sensitivity 3  %12.4le %12.4le %12.4le "%(uS[2][0], uS[2][1], uS[2][2])

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

	nje = cvodes.CVDenseGetNumJacEvals(cvode_mem)
	nfeLS = cvodes.CVDenseGetNumRhsEvals(cvode_mem)

	print "\nFinal Statistics\n"
	print "nst     = %5ld\n"%(nst)
	print "nfe     = %5ld"%(nfe)
	print "netf    = %5ld    nsetups  = %5ld"%(netf, nsetups)
	print "nni     = %5ld    ncfn     = %5ld\n"%(nni, ncfn)

	if sensi:
		print "nfSe    = %5ld    nfeS     = %5ld"%(nfSe, nfeS)
		print "netfs   = %5ld    nsetupsS = %5ld"%(netfS, nsetupsS)
		print "nniS    = %5ld    ncfnS    = %5ld\n"%(nniS, ncfnS)
		print "nje    = %5ld    nfeLS     = %5ld"%(nje, nfeLS)

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

data = UserData()
data.p[0] = 0.04
data.p[1] = 1.0e4
data.p[2] = 3.0e7

y = cvodes.NVector([Y1, Y2, Y3])

cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeMalloc(cvode_mem, f, T0, y, cvodes.CV_WF, 0.0, None)
cvodes.CVodeSetEwtFn(cvode_mem, ewt, None)
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVDense(cvode_mem, NEQ)
cvodes.CVDenseSetJacFn(cvode_mem, Jac, ctypes.pointer(data))

print "\n3-species chemical kinetics problem"

if sensi:
	yStmp = [cvodes.NVector([0,0,0]), cvodes.NVector([0,0,0]), cvodes.NVector([0,0,0])]
	yS = nvecserial.NVectorArray([([0]*NEQ)]*NS)

	cvodes.CVodeSensMalloc(cvode_mem, NS, sensi_meth, yS)
	cvodes.CVodeSetSensRhs1Fn(cvode_mem, fS, ctypes.pointer(data))
	cvodes.CVodeSetSensErrCon(cvode_mem, err_con)
	cvodes.CVodeSetSensParams(cvode_mem, None, data.p, None)

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

print
print "======================================================================="
print "     T     Q       H      NST           y1           y2           y3    "
print "======================================================================="

iout = 1
tout = T1
while iout <= NOUT:
	cvodes.CVode(cvode_mem, tout, y, ctypes.byref(t), cvodes.CV_NORMAL)
	PrintOutput(cvode_mem, t, y)

	if sensi:
		cvodes.CVodeGetSens(cvode_mem, t, yS)
		PrintOutputS(yS)
	print "-----------------------------------------------------------------------"

	iout += 1
	tout *= TMULT

PrintFinalStats(cvode_mem, sensi)
