try:
	from pysundials import ida
	from pysundials import nvecserial
except ImportError:
	import ida
	import nvecserial
import ctypes
import math

NOUT = 11
MGRID = 10
NEQ = MGRID*MGRID
BVAL = 0.1

class UserData(ctypes.Structure):
	_fields_ = [
		("mm", ctypes.c_long),
		("dx", ida.realtype),
		("coeff", ida.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def heatres(tres, uu, up, resval, rdata):
	data = ctypes.cast(rdata, PUserData).contents
	resval[:] = uu
	
	for j in range(1, data.mm-1):
		offset = data.mm*j
		for i in range(1, data.mm-1):
			loc = offset + i
			resval[loc] = up[loc] - data.coeff * (uu[loc-1] + uu[loc+1] + uu[loc-data.mm] + uu[loc+data.mm] - 4.0*uu[loc])
	
	return 0

def SetInitialProfile(data, uu, up, id, res):
	mm1 = data.mm - 1
	
	id = [1.0]*NEQ

	for j in range(data.mm):
		yfact = data.dx * j
		offset = data.mm*j
		for i in range(data.mm):
			xfact = data.dx * i
			loc = offset + i
			uu[loc] = 16.0 * xfact * (1.0 - xfact) * yfact * (1.0 - yfact)
	
	up *= 0

	heatres(0, uu, up, res, ctypes.pointer(data))
	
	up[:] = res * -1

	for j in range(data.mm):
		offset = data.mm*j
		for i in range(data.mm):
			loc = offset + i
			if j == 0 or j == mm1 or i == 0 or i == mm1:
				uu[loc] = BVAL
				up[loc] = 0
				id[loc] = 0
	return 0

def PrintOutput(mem, t, uu):
	umax = max(abs(uu))
	
	kused = ida.IDAGetLastOrder(mem)
	nst = ida.IDAGetNumSteps(mem)
	nni = ida.IDAGetNumNonlinSolvIters(mem)
	nre = ida.IDAGetNumResEvals(mem)
	hused = ida.IDAGetLastStep(mem)
	nje = ida.IDABandGetNumJacEvals(mem)
	nreLS = ida.IDABandGetNumResEvals(mem)

	print " %5.2f %13.5le	%d	%3ld	%3ld	%3ld	%4ld	%4ld	%9.2le "%(t, umax, kused, nst, nni, nje, nre, nreLS, hused)

uu = ida.NVector([0]*NEQ)
up = ida.NVector([0]*NEQ)
res = ida.NVector([0]*NEQ)
constraints = ida.NVector([1]*NEQ)
id = ida.NVector([0]*NEQ)

data = UserData()
data.mm = MGRID
data.dx = 1.0/(MGRID - 1.0)
data.coeff = 1.0/( (data.dx) * (data.dx) )

SetInitialProfile(data, uu, up, id, res)

t0 = 0.0
t1 = 0.01
tret = ida.realtype(0)
rtol = 0
atol = ida.realtype(1.0e-3)

mem = ida.IDACreate()
ida.IDASetRdata(mem, ctypes.pointer(data))
ida.IDASetId(mem, id)
ida.IDASetConstraints(mem, constraints)
ida.IDAMalloc(mem, heatres, t0, uu, up, ida.IDA_SS, rtol, atol)

mu = MGRID
ml = MGRID
ida.IDABand(mem, NEQ, mu, ml)

ida.IDACalcIC(mem, ida.IDA_YA_YDP_INIT, t1)

print "idabanx1: Heat equation, serial example problem for IDA"
print "					Discretized heat equation on 2D unit square."
print "					Zero boundary conditions, polynomial initial conditions."
print "					Mesh dimensions: %i x %i"%(MGRID, MGRID),
print "			Total system size: %d\n"%(NEQ)
print "Tolerance parameters:	rtol = %lg	 atol = %lg"%(rtol, atol.value)
print "Constraints set to force all solution components >= 0. "
print "Linear solver: IDABAND, banded direct solver "
print "			 difference quotient Jacobian, half-bandwidths = %i "%(MGRID)
print "IDACalcIC called with input boundary values = %lg "%(BVAL)
print "\n	 Output Summary (umax = max-norm of solution) \n"
print "	time			 umax		 k	nst	nni	nje	 nre	 nreLS		h			"
print " .	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	. "

PrintOutput(mem, t0, uu)

tout = t1
iout = 1
while iout <= NOUT:
	ida.IDASolve(mem, tout, ctypes.byref(tret), uu, up, ida.IDA_NORMAL)

	PrintOutput(mem, tret.value, uu)
	iout += 1
	tout *= 2.0

netf = ida.IDAGetNumErrTestFails(mem)
ncfn = ida.IDAGetNumNonlinSolvConvFails(mem)
print "\n netf = %ld,	 ncfn = %ld "%(netf, ncfn)
