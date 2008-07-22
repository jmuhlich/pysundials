#!/usr/bin/python

import ctypes
try:
	from pysundials import cvode
except ImportError:
	import cvode
import math

P1_ETA = 3.0
P1_T1 = 1.39283880203
P1_DTOUT = 2.214773875
P1_TOL_FACTOR = 1.0e4

P2_MESHX = 5
P2_MESHY = 5
P2_NEQ = P2_MESHX*P2_MESHY
P2_ALPH1 = 1.0
P2_ALPH2 = 1.0
P2_ML = 5
P2_MU = 0
P2_T0 = 0.0
P2_T1 = 0.01
P2_TOUT_MULT = 10.0
P2_TOL_FACTOR = 1.0e3

def f1(t, y, ydot, f_data):
	ydot[0] = y[1]
	ydot[1] = (1 - y[0]**2) * P1_ETA * y[1] - y[0]
	
	return 0

def Jac1(N, J, tn, y, fy, jac_data, tmp1, tmp2, tmp3):
	J[0][1] = 1.0
	J[1][0] = -2.0 * P1_ETA * y[0] * y[1] - 1.0
	J[1][1] = P1_ETA * (1.0 - y[0]**2)
	
	return 0

def f2(t, y, ydot, f_data):
	for j in range(P2_MESHY):
		for i in range(P2_MESHX):
			k = i + j * P2_MESHX
			d = -2.0*y[k]
			if i != 0:
				d += P2_ALPH1 * y[k-1]
			if j != 0:
				d += P2_ALPH2 * y[k-P2_MESHX]
			ydot[k] = d
	
	return 0

def Jac2(N, mu, ml, J, tn, y, fy, jac_data, tmp1, tmp2, tmp3):
	for j in range(P2_MESHY):
		for i in range(P2_MESHX):
			J[i][j] = -2.0
			if i != P2_MESHX-1:
				J[i][j+1] = P2_ALPH1
			if j != P2_MESHY-1:
				J[i+1][j] = P2_ALPH2
	
	return 0

def MaxError(y, t):
	ex = 0
	maxError = 0
	jfact_inv = 1.0
	
	if t.value == 0:
		return 0
	
	if t.value <= 30:
		ex = math.exp(-2.0*t.value)
	
	for j in range(P2_MESHY):
		ifact_inv = 1.0
		for i in range(P2_MESHX):
			k = i + j * P2_MESHX
			yt = t.value**(i+j) * ex * ifact_inv * jfact_inv
			er = abs(y[k] - yt)
			if er > maxError:
				maxError = er
			ifact_inv /= (i+1)
		jfact_inv /= (j+1)
	return maxError

def PrepareNextRun(cvode_mem, lmm, miter, mu, ml):
	print "\n\n-------------------------------------------------------------"
	print "\nLinear Multistep Method :",
	if lmm == cvode.CV_ADAMS:
		print "ADAMS"
	else:
		print "BDF"
	
	print "Iteration               :",
	if miter == "FUNC":
		print "FUNCTIONAL"
	else:
		print "NEWTON"
		print "Linear Solver           :",
		if miter == "DENSE_USER":
			print "Dense, User-Supplied Jacobian"
			cvode.CVDense(cvode_mem, 2)
			cvode.CVDenseSetJacFn(cvode_mem, Jac1, None)
		elif miter == "DENSE_DQ": 
			print "Dense, Difference Quotient Jacobian"
			cvode.CVDenseSetJacFn(cvode_mem, None, None)
		elif miter == "DIAG": 
			print "Diagonal Jacobian"
			cvode.CVDiag(cvode_mem)
		elif miter == "BAND_USER":
			print "Band, User-Supplied Jacobian"
			cvode.CVBand(cvode_mem, P2_NEQ, mu, ml)
			cvode.CVBandSetJacFn(cvode_mem, Jac2, None)
		elif miter == "BAND_DQ":
			print "Band, Difference Quotient Jacobian"
			cvode.CVBandSetJacFn(cvode_mem, None, None)

def PrintErrOutput(tol_factor):
	print "\n\n Error exceeds %g * tolerance \n"%(tol_factor)

def PrintFinalStats(cvode_mem, miter, ero):
	(lenrw, leniw) = cvode.CVodeGetWorkSpace(cvode_mem)
	nst = cvode.CVodeGetNumSteps(cvode_mem)
	nfe = cvode.CVodeGetNumRhsEvals(cvode_mem)
	nsetups = cvode.CVodeGetNumLinSolvSetups(cvode_mem)
	netf = cvode.CVodeGetNumErrTestFails(cvode_mem)
	nni = cvode.CVodeGetNumNonlinSolvIters(cvode_mem)
	ncfn = cvode.CVodeGetNumNonlinSolvConvFails(cvode_mem)
	
	print "\n Final statistics for this run:\n"
	print " CVode real workspace length              = %4i "%(lenrw.value)
	print " CVode integer workspace length           = %4i "%(leniw.value)
	print " Number of steps                          = %4i "%(nst)
	print " Number of f-s                            = %4i "%(nfe)
	print " Number of setups                         = %4i "%(nsetups)
	print " Number of nonlinear iterations           = %4i "%(nni)
	print " Number of nonlinear convergence failures = %4i "%(ncfn)
	print " Number of error test failures            = %4i \n"%(netf)
	
	if miter != "FUNC":
		if miter == "DENSE_USER" or miter == "DENSE_DQ":
			nje = cvode.CVDenseGetNumJacEvals(cvode_mem)
			nfeLS = cvode.CVDenseGetNumRhsEvals(cvode_mem)
			(lenrwLS, leniwLS) = cvode.CVDenseGetWorkSpace(cvode_mem)
		elif miter == "BAND_USER" or miter == "BAND_DQ":
			nje = cvode.CVBandGetNumJacEvals(cvode_mem)
			nfeLS = cvode.CVBandGetNumRhsEvals(cvode_mem)
			(lenrwLS, leniwLS) = cvode.CVBandGetWorkSpace(cvode_mem)
		elif miter == "DIAG":
			nje = nsetups
			nfeLS = cvode.CVDiagGetNumRhsEvals(cvode_mem)
			(lenrwLS, leniwLS) = cvode.CVDiagGetWorkSpace(cvode_mem)
		print " Linear solver real workspace length      = %4i "%(lenrwLS)
		print " Linear solver integer workspace length   = %4i "%(leniwLS)
		print " Number of Jacobian evaluations           = %4i  "%(nje)
		print " Number of f evals. in linear solver      = %4i \n"%(nfeLS)
	
	print " Error overrun = %.3f "%(ero)

def Problem1():
	nerr = 0
	reltol = cvode.realtype(0)
	abstol = cvode.realtype(1.0e-6)
	t = cvode.realtype(0)

	y = cvode.NVector([0, 0])
	print "Demonstration program for CVODE package - direct linear solvers\n\n"
	print "Problem 1: Van der Pol oscillator"
	print " xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0"
	print " neq = %i,  itol = %s,  reltol = %.2g,  abstol = %.2g"%(2, "CV_SS", 0, 1.0e-6),
	
	cvode_mem = cvode.CVodeCreate(cvode.CV_ADAMS, cvode.CV_FUNCTIONAL)
	
	for miter in ["FUNC", "DENSE_USER", "DENSE_DQ", "DIAG"]:
		ero = 0
		y[0] = 2.0
		y[1] = 0
	
		if miter == "FUNC": 
			cvode.CVodeMalloc(cvode_mem, f1, 0, y, cvode.CV_SS, reltol, abstol)
		else:
			cvode.CVodeSetIterType(cvode_mem, cvode.CV_NEWTON)
			cvode.CVodeReInit(cvode_mem, f1, 0, y, cvode.CV_SS, reltol, abstol)
		
		PrepareNextRun(cvode_mem, cvode.CV_ADAMS, miter, 0, 0)
	
		print "\n     t           x              xdot         qu     hu "
	
		iout = 1
		tout = P1_T1
		while iout <= 4:
			flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
			qu = cvode.CVodeGetLastOrder(cvode_mem)
			hu = cvode.CVodeGetLastStep(cvode_mem)
			print "%10.6g  %14.5e %14.5e   %2i    %6.4e"%(t.value, y[0], y[1], qu, hu)
			if flag != cvode.CV_SUCCESS:
				nerr += 1
				break
			if iout%2 == 0:
				er = abs(y[0])/abstol.value
				if er > ero:
					ero = er
				if er > P1_TOL_FACTOR:
					nerr += 1
					PrintErrOutput(P1_TOL_FACTOR)
			iout += 1
			tout += P1_DTOUT
		
		PrintFinalStats(cvode_mem, miter, ero)
	
	cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_FUNCTIONAL)
	
	for miter in ["FUNC", "DENSE_USER", "DENSE_DQ", "DIAG"]:
		ero = 0
		y[0] = 2.0
		y[1] = 0
	
		if miter == "FUNC": 
			cvode.CVodeMalloc(cvode_mem, f1, 0, y, cvode.CV_SS, reltol, abstol)
		else:
			cvode.CVodeSetIterType(cvode_mem, cvode.CV_NEWTON)
			cvode.CVodeReInit(cvode_mem, f1, 0, y, cvode.CV_SS, reltol, abstol)
			
		PrepareNextRun(cvode_mem, cvode.CV_BDF, miter, 0, 0)
		
		print "\n     t           x              xdot         qu     hu "
			
		iout = 1
		tout = P1_T1
		while iout <= 4:
			flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
			qu = cvode.CVodeGetLastOrder(cvode_mem)
			hu = cvode.CVodeGetLastStep(cvode_mem)
			print "%10.6g  %14.5e %14.5e   %2i    %6.4e"%(t.value, y[0], y[1], qu, hu)
			if flag != cvode.CV_SUCCESS:
				nerr += 1
				break
			if iout%2 == 0:
				er = abs(y[0])/abstol.value
				if er > ero:
					ero = er
				if er > P1_TOL_FACTOR:
					nerr += 1
					PrintErrOutput(P1_TOL_FACTOR)
			iout += 1
			tout += P1_DTOUT
		
		PrintFinalStats(cvode_mem, miter, ero)
	
	return nerr

def Problem2():
	nerr = 0
	reltol = cvode.realtype(0)
	abstol = cvode.realtype(1.0e-6)
	t = cvode.realtype(0)

	print "\n\n-------------------------------------------------------------"
	print "-------------------------------------------------------------"
	print "\nProblem 2: ydot = A * y, where A is a banded lower"
	print "triangular matrix derived from 2-D advection PDE\n"
	print " neq = %i, ml = %i, mu = %i"%(P2_NEQ, P2_ML, P2_MU)
	print " itol = %s, reltol = %.2g, abstol = %.2g"%("CV_SS", 0, 1.0e-6)
	print "      t        max.err      qu     hu "
	
	cvode_mem = cvode.CVodeCreate(cvode.CV_ADAMS, cvode.CV_FUNCTIONAL)
	
	for miter in ["FUNC", "DIAG", "BAND_USER", "BAND_DQ"]:
		ero = 0
		y = cvode.NVector([0]*P2_NEQ)
		y[0] = 1.0
		
		if miter == "FUNC":
			cvode.CVodeMalloc(cvode_mem, f2, P2_T0, y, cvode.CV_SS, reltol, abstol)
		else:
			cvode.CVodeSetIterType(cvode_mem, cvode.CV_NEWTON)
			cvode.CVodeReInit(cvode_mem, f2, P2_T0, y, cvode.CV_SS, reltol, abstol)
		
		PrepareNextRun(cvode_mem, cvode.CV_ADAMS, miter, P2_MU, P2_ML)
		
		print "\n      t        max.err      qu     hu "
		
		iout = 1
		tout = P2_T1
		while iout <= 5:
			flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
			erm = MaxError(y, t)
			qu = cvode.CVodeGetLastOrder(cvode_mem)
			hu = cvode.CVodeGetLastStep(cvode_mem)
			print "%10.3F  %12.4e   %2i   %12.4e"%(t.value, erm, qu, hu)
			if flag != cvode.CV_SUCCESS:
				nerr += 1
				break

			er = erm / abstol.value
			if er > ero:
				ero = er
			if er > P2_TOL_FACTOR:
				nerr += 1
				PrintErrOutput(P2_TOL_FACTOR)
			iout += 1
			tout *= P2_TOUT_MULT
		
		PrintFinalStats(cvode_mem, miter, ero)
	
	cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_FUNCTIONAL)
	
	for miter in ["FUNC", "DIAG", "BAND_USER", "BAND_DQ"]:
		ero = 0
		y[0] = 1.0
		
		if miter == "FUNC":
			cvode.CVodeMalloc(cvode_mem, f2, P2_T0, y, cvode.CV_SS, reltol, abstol)
		else:
			cvode.CVodeSetIterType(cvode_mem, cvode.CV_NEWTON)
			cvode.CVodeReInit(cvode_mem, f2, P2_T0, y, cvode.CV_SS, reltol, abstol)
		
		PrepareNextRun(cvode_mem, cvode.CV_BDF, miter, P2_MU, P2_ML)
		
		print "\n      t        max.err      qu     hu "
		
		iout = 1
		tout = P2_T1
		while iout <= 5:
			flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
			erm = MaxError(y, t)
			qu = cvode.CVodeGetLastOrder(cvode_mem)
			hu = cvode.CVodeGetLastStep(cvode_mem)
			print "%10.3F  %12.4e   %2i   %12.4e"%(t.value, erm, qu, hu)
			if flag != cvode.CV_SUCCESS:
				nerr += 1
				break

			er = erm / abstol.value
			if er > ero:
				ero = er
			if er > P2_TOL_FACTOR:
				nerr += 1
				PrintErrOutput(P2_TOL_FACTOR)
			iout += 1
			tout *= P2_TOUT_MULT
		
		PrintFinalStats(cvode_mem, miter, ero)
	
	return nerr

nerr = Problem1()
nerr += Problem2()
print "\n\n-------------------------------------------------------------"
print "-------------------------------------------------------------"
print "\n Number of errors encountered = %d "%(nerr)
