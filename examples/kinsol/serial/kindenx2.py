try:
	from pysundials import kinsol
	from pysundials import nvecserial
except ImportError:
	import kinsol
	import nvecserial
import ctypes
import math

NVAR = 8
NEQ = 3*NVAR
FTOL = 1.e-5
STOL = 1.e-5

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

def func(y, f, f_data):
	x1 = y[0]; l1 = y[ 8]; u1 = y[16]
	x2 = y[1]; l2 = y[ 9]; u2 = y[17]
	x3 = y[2]; l3 = y[10]; u3 = y[18]
	x4 = y[3]; l4 = y[11]; u4 = y[19]
	x5 = y[4]; l5 = y[12]; u5 = y[20]
	x6 = y[5]; l6 = y[13]; u6 = y[21]
	x7 = y[6]; l7 = y[14]; u7 = y[22]
	x8 = y[7]; l8 = y[15]; u8 = y[23]
	
	eq1 = - 0.1238*x1 + x7 - 0.001637*x2 - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571
	eq2 = 0.2638*x1 - x7 - 0.07745*x2 - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022
	eq3 = 0.3578*x1 + 0.004731*x2 + x6*x8
	eq4 = - 0.7623*x1 + 0.2238*x2 + 0.3461
	eq5 = x1*x1 + x2*x2 - 1
	eq6 = x3*x3 + x4*x4 - 1
	eq7 = x5*x5 + x6*x6 - 1
	eq8 = x7*x7 + x8*x8 - 1
	
	lb1 = l1 - 1.0 - x1
	lb2 = l2 - 1.0 - x2
	lb3 = l3 - 1.0 - x3
	lb4 = l4 - 1.0 - x4
	lb5 = l5 - 1.0 - x5
	lb6 = l6 - 1.0 - x6
	lb7 = l7 - 1.0 - x7
	lb8 = l8 - 1.0 - x8
	
	ub1 = u1 - 1.0 + x1
	ub2 = u2 - 1.0 + x2
	ub3 = u3 - 1.0 + x3
	ub4 = u4 - 1.0 + x4
	ub5 = u5 - 1.0 + x5
	ub6 = u6 - 1.0 + x6
	ub7 = u7 - 1.0 + x7
	ub8 = u8 - 1.0 + x8
	
	f[0] = eq1; f[ 8] = lb1; f[16] = ub1
	f[1] = eq2; f[ 9] = lb2; f[17] = ub2
	f[2] = eq3; f[10] = lb3; f[18] = ub3
	f[3] = eq4; f[11] = lb4; f[19] = ub4
	f[4] = eq5; f[12] = lb5; f[20] = ub5
	f[5] = eq6; f[13] = lb6; f[21] = ub6
	f[6] = eq7; f[14] = lb7; f[22] = ub7
	f[7] = eq8; f[15] = lb8; f[23] = ub8
	
	return 0

def jac(N, J, y, f, jac_data, tmp1, tmp2):
	x1 = y[0]
	x2 = y[1]
	x3 = y[2]
	x4 = y[3]
	x5 = y[4]
	x6 = y[5]
	x7 = y[6]
	x8 = y[7]
	
	J[0][0] = - 0.1238 + 0.004731*x3
	J[0][1] = - 0.001637 - 0.3578*x3
	J[0][2] = 0.004731*x1 - 0.3578*x2
	J[0][3] = - 0.9338
	J[0][6] = 1.0
	
	J[1][0] = 0.2638 + 0.2238*x3
	J[1][1] = - 0.07745 + 0.7623*x3
	J[1][2] = 0.2238*x1 + 0.7623*x2
	J[1][3] = - 0.6734
	J[1][6] = -1.0
	
	J[2][0] = 0.3578
	J[2][1] = 0.004731
	J[2][5] = x8
	J[2][7] = x6
	
	J[3][0] = - 0.7623
	J[3][1] = 0.2238
	
	J[4][0] = 2.0*x1
	J[4][1] = 2.0*x2
	
	J[5][2] = 2.0*x3
	J[5][3] = 2.0*x4
	
	J[6][4] = 2.0*x5
	J[6][5] = 2.0*x6
	
	J[7][6] = 2.0*x7
	J[7][7] = 2.0*x8
	
	for i in range(9):
		J[8+i][i]   = -1.0
		J[8+i][8+i] =  1.0
		
	for i in range(8):
		J[16+i][i]    = 1.0
		J[16+i][16+i] = 1.0
	
	return 0

def PrintOutput(y):
	print "     l=x+1          x         u=1-x"
	print "   ----------------------------------"
	
	for i in range(NVAR):
		print " %10.6lg   %10.6lg   %10.6lg"%(y[i+NVAR], y[i], y[i+2*NVAR])

def PrintFinalStats(kmem):
	nni = kinsol.KINGetNumNonlinSolvIters(kmem)
	nfe = kinsol.KINGetNumFuncEvals(kmem)
	
	nje = kinsol.KINDenseGetNumJacEvals(kmem)
	nfeD = kinsol.KINDenseGetNumFuncEvals(kmem)
	
	print "\nFinal Statistics.. "
	print "nni    = %5ld    nfe   = %5ld "%(nni, nfe)
	print "nje    = %5ld    nfeD  = %5ld "%(nje, nfeD)

print "\nRobot Kinematics Example"
print "8 variables; -1 <= x_i <= 1"
print "KINSOL problem size: 8 + 2*8 = 24 \n"

y = kinsol.NVector([math.sqrt(2)/2.0]*NVAR+[1]*(NVAR*2))
scale = kinsol.NVector([1]*NEQ)
constraints = kinsol.NVector([0]*NVAR+[1]*(NVAR*2))
mset = 1

kmem = kinsol.KINCreate()
kinsol.KINMalloc(kmem, func, y)
kinsol.KINSetConstraints(kmem, constraints)
kinsol.KINSetFuncNormTol(kmem, FTOL)
kinsol.KINSetScaledStepTol(kmem, STOL)
kinsol.KINDense(kmem, NEQ)
kinsol.KINDenseSetJacFn(kmem, jac, None)
kinsol.KINSetMaxSetupCalls(kmem, mset)

print "Initial guess:"
PrintOutput(y);

kinsol.KINSol(kmem, y, kinsol.KIN_LINESEARCH, scale, scale)

print "\nComputed solution:"
PrintOutput(y)

PrintFinalStats(kmem)
