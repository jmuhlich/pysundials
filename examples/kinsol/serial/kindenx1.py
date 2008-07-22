try:
	from pysundials import kinsol
except ImportError:
	import kinsol
import math
import ctypes

class UserData (ctypes.Structure):
	_fields_ = [
		("lb", kinsol.realtype*2),
		("ub", kinsol.realtype*2)
	]
PUserData = ctypes.POINTER(UserData)

def func(u, f, f_data):
	data = ctypes.cast(f_data, PUserData).contents
	f[0] = 0.5 * math.sin(u[0]*u[1]) - 0.25 * u[1] / math.pi - 0.5 * u[0]
	f[1] = (1.0 - 0.25/math.pi)*(math.exp(2.0*u[0])-math.e) + math.e*u[1]/math.pi - 2.0*math.e*u[0]
	f[2] = u[2] - u[0] + data.lb[0]
	f[3] = u[3] - u[0] + data.ub[0]
	f[4] = u[4] - u[1] + data.lb[1]
	f[5] = u[5] - u[1] + data.ub[1]
	return 0

def SetInitialGuess1(u, data):
	u[0] = data.lb[0]
	u[1] = data.lb[1]
	u[2] = data.lb[0] - data.lb[0]
	u[3] = data.lb[0] - data.ub[0]
	u[4] = data.lb[1] - data.lb[1]
	u[5] = data.lb[1] - data.ub[1]

def SetInitialGuess2(u, data):
	x1 = 0.5 * (data.lb[0] + data.ub[0])
	x2 = 0.5 * (data.lb[1] + data.ub[1])

	u[0] = x1
	u[1] = x2
	u[2] = x1 - data.lb[0]
	u[3] = x1 - data.ub[0]
	u[4] = x2 - data.lb[1]
	u[5] = x2 - data.ub[1]

def PrintFinalStats(kmem):
	nni = kinsol.KINGetNumNonlinSolvIters(kmem)
	nfe = kinsol.KINGetNumFuncEvals(kmem)
	nje = kinsol.KINDenseGetNumJacEvals(kmem)
	nfeD = kinsol.KINDenseGetNumFuncEvals(kmem)
	print "Final Statistics:"
	print "  nni = %5ld    nfe  = %5ld "%(nni, nfe)
	print "  nje = %5ld    nfeD = %5ld "%(nje, nfeD)

def SolveIt(kmem, u, s, glstr, mset):
	print
  	
	if mset == 1:
		print "Exact Newton",
	else:
		print "Modified Newton",

	if glstr == kinsol.KIN_NONE:
		print
	else:
		print "with line search"

	kinsol.KINSetMaxSetupCalls(kmem, mset)
	kinsol.KINSol(kmem, u, glstr, s, s)
	
	print "Solution:\n  [x1,x2] = %r"%(u)
	PrintFinalStats(kmem)
	
	return 0

data = UserData()
data.lb[0] = 0.25
data.lb[1] = 1.5
data.ub[0] = 1.0
data.ub[1] = 2*math.pi

u1 = kinsol.NVector([0.0]*6)
u2 = kinsol.NVector([0.0]*6)
u = kinsol.NVector([0.0]*6)
s = kinsol.NVector([1.0]*6)

SetInitialGuess1(u1,data)
SetInitialGuess2(u2,data)

c = kinsol.NVector([0.0, 0.0, 1.0, -1.0, 1.0, -1.0])
  
fnormtol = 1.0e-5
scsteptol = 1.0e-5

kmem = kinsol.KINCreate()

kinsol.KINSetFdata(kmem, ctypes.pointer(data))
kinsol.KINSetConstraints(kmem, c)
kinsol.KINSetFuncNormTol(kmem, fnormtol)
kinsol.KINSetScaledStepTol(kmem, scsteptol)
kinsol.KINMalloc(kmem, func, u)
kinsol.KINDense(kmem, 6)

print "\nFerraris and Tronconi test problem"
print "Tolerance parameters:"
print "  fnormtol  = %10.6lg\n  scsteptol = %10.6lg"%(fnormtol, scsteptol)

print "\n------------------------------------------"
print "\nInitial guess on lower bounds"
print "  [x1,x2] = %r"%(u1)

#u = 1.0*u1
for i in range(len(u1)):
	u[i] = u1[i]
glstr = kinsol.KIN_NONE
mset = 1
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u1)):
	u[i] = u1[i]
glstr = kinsol.KIN_LINESEARCH
mset = 1
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u1)):
	u[i] = u1[i]
glstr = kinsol.KIN_NONE
mset = 0
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u1)):
	u[i] = u1[i]
glstr = kinsol.KIN_LINESEARCH
mset = 0
SolveIt(kmem, u, s, glstr, mset)

print "------------------------------------------"
print "Initial guess in middle of feasible region"
print "  [x1,x2] = %r"%(u2)

#u = 1.0*u1
for i in range(len(u2)):
	u[i] = u2[i]
glstr = kinsol.KIN_NONE
mset = 1
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u2)):
	u[i] = u2[i]
glstr = kinsol.KIN_LINESEARCH
mset = 1
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u2)):
	u[i] = u2[i]
glstr = kinsol.KIN_NONE
mset = 0
SolveIt(kmem, u, s, glstr, mset)

#u = 1.0*u1
for i in range(len(u2)):
	u[i] = u2[i]
glstr = kinsol.KIN_LINESEARCH
mset = 0
SolveIt(kmem, u, s, glstr, mset)
