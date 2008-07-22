try:
	from pysundials import ida
	from pysundials import nvecserial
except ImportError:
	import ida
	import nvecserial
import math
import ctypes

NOUT = 11
MGRID = 10
NEQ = MGRID*MGRID

class UserData(ctypes.Structure):
	_fields_ = [
		('mm', ctypes.c_long),
		('dx', ida.realtype),
		('coeff', ida.realtype),
		('pp', nvecserial.PVector)
	]
PUserData = ctypes.POINTER(UserData)

def resHeat(tt, uu, up, rr, res_data):
	data = ctypes.cast(res_data, PUserData).contents
	
	rr[:] = uu
	
	for j in range(1, MGRID-1):
		offset = data.mm*j
		for i in range(1, data.mm-1):
			loc = offset + i
			dif1 = uu[loc-1]  + uu[loc+1]  - 2.0 * uu[loc]
			dif2 = uu[loc-data.mm] + uu[loc+data.mm] - 2.0 * uu[loc]
			rr[loc] = up[loc] - data.coeff * ( dif1 + dif2 )
	
	return 0

def PsetupHeat(tt, uu, up, rr, c_j, prec_data, tmp1, tmp2, tmp3):
	data = ctypes.cast(prec_data, PUserData).contents
	
	ppv = ida.NVector(data.pp)
	
	ppv[:] = [1]*len(ppv)
	
	pelinv = 1.0/(c_j + 4.0*data.coeff)
	
	for j in range(1, data.mm-1):
		offset = data.mm * j;
		for i in range(1,data.mm-1):
			loc = offset + i
			ppv[loc] = pelinv
	
	return 0

def PsolveHeat(tt, uu, up, rr, rvec, zvec, c_j, delta, prec_data, tmp):
	data = ctypes.cast(prec_data, PUserData).contents
	zvec[:] = ida.NVector(data.pp)*rvec
	return 0 

def SetInitialProfile(data, uu, up, res):
	mm1 = data.mm - 1
	for j in range(data.mm):
		yfact = data.dx * j
		offset = data.mm*j
		for i in range(data.mm):
			xfact = data.dx * i
			loc = offset + i
			uu[loc] = 16.0 * xfact * (1.0 - xfact) * yfact * (1.0 - yfact)
	
	up[:] = [0]*len(up)
	
	resHeat(0.0, uu, up, res, ctypes.pointer(data))
	
	up[:] = -res
	
	for j in range(data.mm):
		offset = data.mm*j
		for i in range(data.mm):
			loc = offset + i
			if j == 0 or j == mm1 or i == 0 or i == mm1:
				up[loc] = 0.0
	
	return 0

def PrintOutput(mem, t, uu):
	#umax = N_VMaxNorm(uu);
	umax = max(abs(uu))
	
	kused = ida.IDAGetLastOrder(mem)
	nst = ida.IDAGetNumSteps(mem)
	nni = ida.IDAGetNumNonlinSolvIters(mem)
	nre = ida.IDAGetNumResEvals(mem)
	hused = ida.IDAGetLastStep(mem)
	nje = ida.IDASpilsGetNumJtimesEvals(mem)
	nli = ida.IDASpilsGetNumLinIters(mem)
	nreLS = ida.IDASpilsGetNumResEvals(mem)
	npe = ida.IDASpilsGetNumPrecEvals(mem)
	nps = ida.IDASpilsGetNumPrecSolves(mem)
	
	print " %5.2f %13.5le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2le  %3ld %3ld"%(t.value, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps)

tret = ida.realtype(0)
uu = ida.NVector([0]*NEQ)
up = ida.NVector([0]*NEQ)
res = ida.NVector([0]*NEQ)

constraints = ida.NVector([1]*NEQ)

data = UserData()

data.mm  = MGRID
data.dx = 1.0/(MGRID-1.0)
data.coeff = 1.0/(data.dx * data.dx)
params = ida.NVector([0]*NEQ)
data.pp = params.data

SetInitialProfile(data, uu, up, res)

t0 = 0.0
t1 = 0.01
rtol = 0.0
atol = 1.0e-3

mem = ida.IDACreate()
ida.IDASetRdata(mem, ctypes.pointer(data))
ida.IDASetConstraints(mem, constraints)
ida.IDAMalloc(mem, resHeat, t0, uu, up, ida.IDA_SS, rtol, atol)
ida.IDASpgmr(mem, 0)
ida.IDASpilsSetPreconditioner(mem, PsetupHeat, PsolveHeat, ctypes.pointer(data))

print "\nidakryx: Heat equation, serial example problem for IDA "
print "         Discretized heat equation on 2D unit square. "
print "         Zero boundary conditions,",
print " polynomial initial conditions."
print "         Mesh dimensions: %d x %d"%(MGRID, MGRID)
print "        Total system size: %d\n"%(NEQ)
print "Tolerance parameters:  rtol = %lg   atol = %lg"%(rtol, atol)
print "Constraints set to force all solution components >= 0. "
print "Linear solver: IDASPGMR, preconditioner using diagonal elements. "

print "\n\nCase 1: gsytpe = MODIFIED_GS"
print "\n   Output Summary (umax = max-norm of solution) \n"
print "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps"
print "----------------------------------------------------------------------"

tout = t1
iout = 1
while iout <= NOUT:
	ida.IDASolve(mem, tout, ctypes.byref(tret), uu, up, ida.IDA_NORMAL)
	PrintOutput(mem, tret, uu)
	iout += 1
	tout *= 2

netf = ida.IDAGetNumErrTestFails(mem)
ncfn = ida.IDAGetNumNonlinSolvConvFails(mem)
ncfl = ida.IDASpilsGetNumConvFails(mem)

print  "\nError test failures            = %ld"%(netf)
print  "Nonlinear convergence failures = %ld"%(ncfn)
print  "Linear convergence failures    = %ld"%(ncfl)

SetInitialProfile(data, uu, up, res)

ida.IDAReInit(mem, resHeat, t0, uu, up, ida.IDA_SS, rtol, atol)
ida.IDASpilsSetGSType(mem, ida.CLASSICAL_GS)

print "\n\nCase 2: gstype = CLASSICAL_GS"
print "\n   Output Summary (umax = max-norm of solution) \n"
print "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps"
print "----------------------------------------------------------------------"

tout = t1
iout = 1
while iout <= NOUT:
	ida.IDASolve(mem, tout, ctypes.byref(tret), uu, up, ida.IDA_NORMAL)
	PrintOutput(mem, tret, uu)
	iout += 1
	tout *= 2

netf = ida.IDAGetNumErrTestFails(mem)
ncfn = ida.IDAGetNumNonlinSolvConvFails(mem)
ncfl = ida.IDASpilsGetNumConvFails(mem)

print  "\nError test failures            = %ld"%(netf)
print  "Nonlinear convergence failures = %ld"%(ncfn)
print  "Linear convergence failures    = %ld"%(ncfl)
