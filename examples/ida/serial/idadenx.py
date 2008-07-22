try:
	from pysundials import ida
except ImportError:
	import ida
import ctypes

def resrob(tres, yy, yp, rr, rdata):
	rr[0]  = -0.04*yy[0] + 1.0e4*yy[1]*yy[2]
	rr[1]  = -rr[0] - 3.0e7*yy[1]*yy[1] - yp[1]
	rr[0] -=  yp[0]
	rr[2]  =  yy[0] + yy[1] + yy[2] - 1.0
	return 0

def grob(t, yy, yp, gout, g_data):
	gout[0] = yy[0] - 0.0001
	gout[1] = yy[2] - 0.01
	return 0

def jacrob(Neq, tt, yy, yp, resvec, cj, jdata, JJ, tempv1, tempv2, tempv3):
	JJ[0][0] = -0.04 - cj
	JJ[1][0] = 0.04
	JJ[2][0] = 1.0
	JJ[0][1] = 1.0e4*yy[2]
	JJ[1][1] = -1.0e4*yy[2] - 6.0e7*yy[1] - cj
	JJ[2][1] = 1.0
	JJ[0][2] = 1.0e4*yy[1]
	JJ[1][2] = -1.0e4*yy[1]
	JJ[2][2] = 1.0
	
	return 0

def PrintOutput(mem, t, y):
	kused = ida.IDAGetLastOrder(mem)
	nst = ida.IDAGetNumSteps(mem)
	hused = ida.IDAGetLastStep(mem)
	print "%10.4le %12.4le %12.4le %12.4le | %3ld  %1d %12.4le\n"%(t.value, y[0], y[1], y[2], nst, kused, hused)

def PrintRootInfo(root_f1, root_f2):
	print "    rootsfound[] = %3d %3d"%(root_f1, root_f2)
	return

yy = ida.NVector([1.0, 0.0, 0.0])
yp = ida.NVector([-0.04, 0.04, 0.0])
rtol = 1.0e-4
avtol = ida.NVector([1.0e-8, 1.0e-14, 1.0e-6])

t0 = 0.0
tout1 = 0.4
tret = ida.realtype(0.0)

print"\nidadenx: Robertson kinetics DAE serial example problem for IDA"
print"         Three equation chemical kinetics problem.\n"
print"Linear solver: IDADENSE, with user-supplied Jacobian."
print"Tolerance parameters:  rtol = %lg   atol = %lg %lg %lg \n"%(rtol, avtol[0], avtol[1], avtol[2])
print"Initial conditions y0 = (%lg %lg %lg)"%(yy[0], yy[1], yy[2])
print"Constraints and id not used.\n"
print"-----------------------------------------------------------------------"
print"  t             y1           y2           y3      | nst  k      h"
print"-----------------------------------------------------------------------"


mem = ida.IDACreate()
ida.IDAMalloc(mem, resrob, t0, yy, yp, ida.IDA_SV, rtol, avtol)
del avtol
ida.IDARootInit(mem, 2, grob, None)
ida.IDADense(mem, 3)
ida.IDADenseSetJacFn(mem, jacrob, None)

iout = 0
tout = tout1

while True:
	retval = ida.IDASolve(mem, tout, ctypes.byref(tret), yy, yp, ida.IDA_NORMAL)
	PrintOutput(mem,tret,yy)
	
	if retval == ida.IDA_ROOT_RETURN:
	  rootsfound = ida.IDAGetRootInfo(mem, 2)
	  PrintRootInfo(rootsfound[0],rootsfound[1])
	elif retval == ida.IDA_SUCCESS:
	  iout += 1
	  tout *= 10.0
	else:
		print retval
		break

	
	
	if iout == 12:
		break
  
	print "\nFinal Run Statistics: \n"
	print "Number of steps                    = %ld"%ida.IDAGetNumSteps(mem)
	print "Number of residual evaluations     = %ld"%(ida.IDAGetNumResEvals(mem)+ida.IDADenseGetNumResEvals(mem))
	print "Number of Jacobian evaluations     = %ld"%ida.IDADenseGetNumJacEvals(mem)
	print "Number of nonlinear iterations     = %ld"%ida.IDAGetNumNonlinSolvIters(mem)
	print "Number of error test failures      = %ld"%ida.IDAGetNumErrTestFails(mem)
	print "Number of nonlinear conv. failures = %ld"%ida.IDAGetNumNonlinSolvConvFails(mem)
	print "Number of root fn. evaluations     = %ld"%ida.IDAGetNumGEvals(mem)
