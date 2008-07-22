#!/usr/bin/python

try:
	from pysundials import cvode
except ImportError:
	import cvode
import ctypes

def f(t, y, ydot, f_data):
	yd1 = ydot[0] = -0.04*y[0] + 1.0e4*y[1]*y[2]
	yd3 = ydot[2] = 3.0e7*y[1]*y[1]
	ydot[1] = -yd1 - yd3
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

y = cvode.NVector([1.0, 0.0, 0.0])
abstol = cvode.NVector([1.0e-8, 1.0e-14, 1.0e-6])
reltol = cvode.realtype(1.0e-4)

cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SV, reltol, abstol)
cvode.CVodeRootInit(cvode_mem, 2, g, None)
cvode.CVDense(cvode_mem, 3)
cvode.CVDenseSetJacFn(cvode_mem, Jac, None)

print " \n3-species kinetics problem\n"

iout = 0
tout = 0.4
t = cvode.realtype(0.0)

while True:
	flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
	print "At t = %-14.4e  y =  %-11.6e    %-11.6e    %-11.6e"%(t.value, y[0], y[1], y[2])
	
	if flag == cvode.CV_ROOT_RETURN:
		rootsfound = cvode.CVodeGetRootInfo(cvode_mem, 2)
		print "    rootsfound[] = %3i %3i"%(rootsfound[0],rootsfound[1])
	
	if flag == cvode.CV_SUCCESS:
		iout += 1;
		tout *= 10.0;

	if iout == 12:
		break

nst = cvode.CVodeGetNumSteps(cvode_mem)
nfe = cvode.CVodeGetNumRhsEvals(cvode_mem)
nsetups = cvode.CVodeGetNumLinSolvSetups(cvode_mem)
netf = cvode.CVodeGetNumErrTestFails(cvode_mem)
nni = cvode.CVodeGetNumNonlinSolvIters(cvode_mem)
ncfn = cvode.CVodeGetNumNonlinSolvConvFails(cvode_mem)
nje = cvode.CVDenseGetNumJacEvals(cvode_mem)
nfeLS = cvode.CVDenseGetNumRhsEvals(cvode_mem)
nge = cvode.CVodeGetNumGEvals(cvode_mem)

print "\nFinal Statistics:"
print "nst = %-6i nfe  = %-6i nsetups = %-6i nfeLS = %-6i nje = %i"%(nst, nfe, nsetups, nfeLS, nje)
print "nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n "%(nni, ncfn, netf, nge)
