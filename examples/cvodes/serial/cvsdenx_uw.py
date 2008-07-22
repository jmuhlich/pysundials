#!/usr/bin/python

try:
	from pysundials import cvodes
except ImportError:
	import cvodes
import ctypes

def f(t, y, ydot, f_data):
	yd1 = ydot[0] = -0.04*y[0] + 1.0e4*y[1]*y[2]
	yd3 = ydot[2] = 3.0e7*y[1]*y[1]
	ydot[1] = -yd1 - yd3
	return 0

def ewt(y, w, e_data):
	atol = [1.0e-8, 1.0e-14, 1.0e-6]
	
	for i in range(3):
		yy = y[i]
		ww = 1.0e-4*abs(yy) + atol[i]
		if ww < 0.0:
			return -1
		w[i] = 1/ww
	
	return 0

def g(t, y, gout, g_data):
	gout[0] = y[0] - 0.0001
	gout[1] = y[2] - 0.01
	return 0

def Jac(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):
	J[0][0] = -0.04
	J[0][1] = 1.0e4*y[2]
	J[0][2] = 1.0e4*y[1]
	J[1][0] = 0.04
	J[1][1] = -1.0e4*y[2]-6.0e7*y[1]
	J[1][2] = -1.0e4*y[1]
	J[2][1] = 6.0e7*y[1]
	return 0

y = cvodes.NVector([1.0, 0.0, 0.0])

cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeMalloc(cvode_mem, f, 0.0, y, cvodes.CV_WF, 0.0, None)
cvodes.CVodeSetEwtFn(cvode_mem, ewt, None)
cvodes.CVodeRootInit(cvode_mem, 2, g, None)
cvodes.CVDense(cvode_mem, 3)
cvodes.CVDenseSetJacFn(cvode_mem, Jac, None)

print " \n3-species kinetics problem\n"

iout = 0
tout = 0.4
t = cvodes.realtype(0.0)

while True:
	flag = cvodes.CVode(cvode_mem, tout, y, ctypes.byref(t), cvodes.CV_NORMAL)
	print "At t = %-14.4e  y =  %-11.6e    %-11.6e    %-11.6e"%(t.value, y[0], y[1], y[2])
	
	if flag == cvodes.CV_ROOT_RETURN:
		rootsfound = cvodes.CVodeGetRootInfo(cvode_mem,2)
		print "    rootsfound[] = %3i %3i"%(rootsfound[0],rootsfound[1])
	
	if flag == cvodes.CV_SUCCESS:
		iout += 1;
		tout *= 10.0;

	if iout == 12:
		break

nst = cvodes.CVodeGetNumSteps(cvode_mem)
nfe = cvodes.CVodeGetNumRhsEvals(cvode_mem)
nsetups = cvodes.CVodeGetNumLinSolvSetups(cvode_mem)
netf = cvodes.CVodeGetNumErrTestFails(cvode_mem)
nni = cvodes.CVodeGetNumNonlinSolvIters(cvode_mem)
ncfn = cvodes.CVodeGetNumNonlinSolvConvFails(cvode_mem)
nje = cvodes.CVDenseGetNumJacEvals(cvode_mem)
nfeLS = cvodes.CVDenseGetNumRhsEvals(cvode_mem)
nge = cvodes.CVodeGetNumGEvals(cvode_mem)

print "\nFinal Statistics:"
print "nst = %-6i nfe  = %-6i nsetups = %-6i nfeLS = %-6i nje = %i"%(nst, nfe, nsetups, nfeLS, nje)
print "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n "%(nni, ncfn, netf, nge)
