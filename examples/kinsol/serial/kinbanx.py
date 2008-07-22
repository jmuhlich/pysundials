try:
	from pysundials import kinsol
	from pysundials import nvecserial
except ImportError:
	import kinsol
	import nvecserial
import ctypes
import math

NX = 31
NY = 31
NEQ = NX*NY
SKIP = 3
FTOL = 1.e-12

#define IJth(vdata,i,j) (vdata[(j-1) + (i-1)*NY])

def func(u, f, f_data):
	dx = 1.0/(NX+1)
	dy = 1.0/(NY+1)
	hdc = 1.0/(dx*dx)
	vdc = 1.0/(dy*dy)
	
	for j in range(NY):
		y = (j+1)*dy
		for i in range(NX):
			x = (i+1)*dx
			uij = u[(j) + (i)*NY]
			if j == 0:
				udn = 0.0
			else:
				udn = u[(j-1) + (i)*NY]
			if j == NY-1:
				uup = 0.0
			else:
				uup = u[(j+1) + (i)*NY]
			if i == 0:
				ult = 0.0
			else:
				ult = u[(j) + (i-1)*NY]
			if i == NX-1:
				urt = 0.0
			else:
				urt = u[(j) + (i+1)*NY]
			hdiff = hdc*(ult - 2.0*uij + urt)
			vdiff = vdc*(uup - 2.0*uij + udn)
			f[(j) + (i)*NY] = hdiff + vdiff + uij - uij*uij*uij + 2.0
	
	return 0

def PrintOutput(u):
	dx = 1.0/(NX+1)
	dy = 1.0/(NY+1)
	
	print "            "
	for i in range(1, NX+1, SKIP):
		x = i*dx
		print "%-8.5lf"%(x),
		
	print  "\n"
	
	for j in range(1, NY+1, SKIP):
		y = j*dy
		print "%-8.5lf   "%(y),
		for i in range(1, NX+1, SKIP):
			print "%-8.5lf"%(u[(j-1) + (i-1)*NY]),
		print

def PrintFinalStats(kmem):
	nni = kinsol.KINGetNumNonlinSolvIters(kmem)
	nfe = kinsol.KINGetNumFuncEvals(kmem)
	nbcfails = kinsol.KINGetNumBetaCondFails(kmem)
	nbacktr = kinsol.KINGetNumBacktrackOps(kmem)
	
	(lenrw, leniw) = kinsol.KINGetWorkSpace(kmem)
	
	nje = kinsol.KINBandGetNumJacEvals(kmem)
	nfeD = kinsol.KINBandGetNumFuncEvals(kmem)
	
	(lenrwB, leniwB) = kinsol.KINBandGetWorkSpace(kmem)
	
	print "\nFinal Statistics.. \n"
	print "nni      = %6ld    nfe     = %6ld "%(nni, nfe)
	print "nbcfails = %6ld    nbacktr = %6ld "%(nbcfails, nbacktr)
	print "nje      = %6ld    nfeB    = %6ld "%(nje, nfeD)
	print "lenrw    = %6ld    leniw   = %6ld "%(lenrw, leniw)
	print "lenrwB   = %6ld    leniwB  = %6ld "%(lenrwB, leniwB)

print "\n2D elliptic PDE on unit square"
print "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0"
print " + homogeneous Dirichlet boundary conditions\n"
print "Solution method: Modified Newton with band linear solver"
print "Problem size: %2ld x %2ld = %4ld"%(NX, NY, NEQ)

y = kinsol.NVector([0]*NEQ)
scale = kinsol.NVector([1]*NEQ)

mset = 100
msubset = 1

kmem = kinsol.KINCreate()
kinsol.KINMalloc(kmem, func, y)
kinsol.KINSetFuncNormTol(kmem, FTOL)
kinsol.KINBand(kmem, NEQ, NY, NY)
kinsol.KINSetMaxSetupCalls(kmem, mset)
kinsol.KINSetMaxSubSetupCalls(kmem, msubset)
kinsol.KINSol(kmem, y, kinsol.KIN_LINESEARCH, scale, scale)
fnorm = kinsol.KINGetFuncNorm(kmem)

print "\nComputed solution (||F|| = %g):\n"%(fnorm)
PrintOutput(y)
PrintFinalStats(kmem)
