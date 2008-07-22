try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math
import sys

XMAX = 2.0
MX = 10
NEQ = MX
ATOL = 1.e-5
T0 = 0.0
T1 = 0.5
DTOUT = 0.5
NOUT = 10

NP = 2
NS = 2

class UserData(ctypes.Structure):
	_fields_ = [
		('p', cvodes.realtype*NP),
		('dx', cvodes.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def f(t, u, udot, f_data):
	data = ctypes.cast(f_data, PUserData).contents
	dx = data.dx
	hordc = data.p[0]/(dx*dx)
	horac = data.p[1]/(2.0*dx)

	for i in range(NEQ):
		ui = u[i]
		if i != 0:
			ult = u[i-1]
		else:
			ult = 0
		if i != NEQ-1:
			urt = u[i+1]
		else:
			urt = 0

		hdiff = hordc*(ult - 2.0*ui + urt)
		hadv = horac*(urt - ult)
		udot[i] = hdiff + hadv
	
	return 0

def WrongArgs(name):
	print "\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]"%(name)
	print "         sensi_meth = sim, stg, or stg1"
	print "         err_con    = t or f"

	sys.exit(0)

def SetIC(u, dx):
	for i in range(NEQ):
		x = (i+1)*dx
		u[i] = x*(XMAX - x)*math.exp(2.0*x)

def PrintOutput(cvode_mem, t, u):
	nst = cvodes.CVodeGetNumSteps(cvode_mem)
	qu = cvodes.CVodeGetLastOrder(cvode_mem)
	hu = cvodes.CVodeGetLastStep(cvode_mem)

	print "%8.3le %2d  %8.3le %5ld"%(t.value, qu, hu, nst)
	print "                                Solution      ",

	print "%12.4le "%(max(abs(u)))

def PrintOutputS(uS):
	print "                                Sensitivity 1 ",
	print "%12.4le "%(max(abs(uS[0])))
	print "                                Sensitivity 2 ",
	print "%12.4le "%(max(abs(uS[1])))

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

	print "\nFinal Statistics\n"
	print "nst     = %5ld\n"%(nst)
	print "nfe     = %5ld"%(nfe)
	print "netf    = %5ld    nsetups  = %5ld"%(netf, nsetups)
	print "nni     = %5ld    ncfn     = %5ld\n"%(nni, ncfn)

	if sensi:
		print "nfSe    = %5ld    nfeS     = %5ld"%(nfSe, nfeS)
		print "netfs   = %5ld    nsetupsS = %5ld"%(netfS, nsetupsS)
		print "nniS    = %5ld    ncfnS    = %5ld"%(nniS, ncfnS)

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
dx = data.dx = XMAX/(float(MX)+1)
data.p[0] = 1.0
data.p[1] = 0.5

u = cvodes.NVector([0]*NEQ)
SetIC(u, dx)

reltol = 0
abstol = ATOL

cvode_mem = cvodes.CVodeCreate(cvodes.CV_ADAMS, cvodes.CV_FUNCTIONAL)
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVodeMalloc(cvode_mem, f, T0, u, cvodes.CV_SS, reltol, abstol)

print "\n1-D advection-diffusion equation, mesh size =%3d"%(MX)

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

print
print "============================================================"
print "     T     Q       H      NST                    Max norm   "
print "============================================================"

iout = 1
tout = T1
while iout <= NOUT:
	cvodes.CVode(cvode_mem, tout, u, ctypes.byref(t), cvodes.CV_NORMAL)
	PrintOutput(cvode_mem, t, u)

	if sensi:
		cvodes.CVodeGetSens(cvode_mem, t, uS)
		PrintOutputS(uS)
	print "------------------------------------------------------------"

	iout += 1
	tout += DTOUT

PrintFinalStats(cvode_mem, sensi)
