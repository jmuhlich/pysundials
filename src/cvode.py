#cvode.py is part of the PySUNDIALS package, and is released under the
#following terms and conditions.

#Copyright (c) 2007, James Dominy, Brett Olivier, Jan Hendrik Hofmeyr, Johann Rohwer
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. Neither the name of the <ORGANIZATION> nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.

"""Python bindings for the CVODE integrator

The following SUNDIALS header files are wrapped:
	cvode.h, cvode_band.h, cvode_bandpre.h, cvode_bbdpre.h, cvode_dense.h,
	cvode_diag.h, cvode_spbcgs.h, cvode_spmgr.h, cvode_spils.h,
	cvode_sptfqmr.h
"""

import ctypes
import sundials_core
import nvecserial

realtype = nvecserial.realtype
NVector = nvecserial.NVector

cvode = sundials_core.loadlib("cvode")

#Linear Multistep method constants
CV_ADAMS = 1
CV_BDF = 2

#Internal time step nonlinear solver
CV_FUNCTIONAL = 1
CV_NEWTON = 2

#Relative and Absolute Tolerance Types
CV_SS = 1
CV_SV = 2
CV_WF = 3

#Internal task job type
CV_NORMAL = 1
CV_ONE_STEP = 2
CV_NORMAL_TSTOP = 3
CV_ONE_STE_TSTOP = 4

#CVode() return flags
CV_SUCCESS = 0
CV_TSTOP_RETURN = 1
CV_ROOT_RETURN = 2

CV_WARNING = 99

CV_TOO_MUCH_WORK = -1
CV_TOO_MUCH_ACC = -2
CV_ERR_FAILURE = -3
CV_CONV_FAILURE = -4

CV_LINIT_FAIL = -5
CV_LSETUP_FAIL = -6
CV_LSOLVE_FAIL = -7
CV_RHSFUNC_FAIL = -8
CV_FIRST_RHSFUNC_ERR = -9
CV_REPTD_RHSFUNC_ERR = -10
CV_UNREC_RHSFUNC_ERR = -11
CV_RTFUNC_FAIL = -12

CV_MEM_FAIL = -20
CV_MEM_NULL = -21
CV_ILL_INPUT = -22
CV_NO_MALLOC = -23
CV_BAD_K = -24
CV_BAD_T = -25
CV_BAD_DKY = -26

#convfail return constants
CV_NO_FAILURES = 0
CV_FAIL_BAD_J = 1
CV_FAIL_OTHER = 2

__Callback = []
__ActualCallback = []

#it's here if I ever need it for something
#class CVodeMemRec(ctypes.Structure):
#	_fields_ = [
#		('cv_uround', realtype),
#		('cv_f', ctypes.c_void_p),
#		('cv_f_data', ctypes.c_void_p),
#		('cv_lmm', ctypes.c_int),
#		('cv_iter', ctypes.c_int),
#		('cv_itol', ctypes.c_int),
#		('cv_reltol', realtype),
#		('cv_Sabstol', realtype),
#		('cv_Vabstol', nvecserial.PVector),
#		('cv_efun', ctypes.c_void_p),
#		('cv_e_data', ctypes.c_void_p),
#		('cv_zn', ctypes.c_void_p),
#		('cv_ewt', nvecserial.PVector),
#		('cv_y', nvecserial.PVector),
#		('cv_acor', nvecserial.PVector),
#		('cv_tempv', nvecserial.PVector),
#		('cv_ftemp', nvecserial.PVector),
#		('cv_tstopset', ctypes.c_int),
#		('cv_istop', ctypes.c_int),
#		('cv_tstop', realtype),
#		('cv_q', ctypes.c_int),
#		('cv_qprime', ctypes.c_int), 
#		('cv_next_q', ctypes.c_int),
#		('cv_qwait', ctypes.c_int),
#		('cv_L', ctypes.c_int),
#		('cv_hin', realtype),
#		('cv_h', realtype),
#		('cv_hprime', realtype), 
#		('cv_next_h', realtype), 
#		('cv_eta', realtype),
#		('cv_hscale', realtype),
#		('cv_tn', realtype),
#		('cv_tretlast', realtype),
#		('cv_tau', ctypes.c_void_p),
#		('cv_tq', ctypes.c_void_p),
#		('cv_l', ctypes.c_void_p),
#		('cv_rl1', realtype),
#		('cv_gamma', realtype),
#		('cv_gammap', realtype),
#		('cv_gamrat', realtype),
#		('cv_crate', realtype),
#		('cv_acnrm', realtype),
#		('cv_nlscoef', realtype),
#		('cv_mnewt', ctypes.c_int),
#		('cv_qmax', ctypes.c_int),
#		('cv_mxstep', ctypes.c_long),
#		('cv_maxcor', ctypes.c_int),
#		('cv_mxhnil', ctypes.c_int),
#		('cv_maxnef', ctypes.c_int),
#		('cv_maxncf', ctypes.c_int),
#		('cv_hmin', realtype),
#		('cv_hmax_inv', realtype),
#		('cv_etamax', realtype),
#		('cv_nst', ctypes.c_long),
#		('cv_nfe', ctypes.c_long),
#		('cv_ncfn', ctypes.c_long),
#		('cv_netf', ctypes.c_long),
#		('cv_nni', ctypes.c_long),
#		('cv_nsetups', ctypes.c_long),
#		('cv_nhnil', ctypes.c_int),
#		('cv_etaqm1', realtype),
#		('cv_etaq', realtype),
#		('cv_etaqp1', realtype),
#		('cv_lrw1', ctypes.c_long), 
#		('cv_liw1', ctypes.c_long), 
#		('cv_lrw', ctypes.c_long),
#		('cv_liw', ctypes.c_long),
#		('cv_linit', ctypes.c_void_p),
#		('cv_lsetup', ctypes.c_void_p), 
#		('cv_lsolve', ctypes.c_void_p),
#		('cv_lfree', ctypes.c_void_p),
#		('cv_lmem', ctypes.c_void_p),           
#		('cv_qu', ctypes.c_int),
#		('cv_nstlp', ctypes.c_long),
#		('cv_h0u', realtype),
#		('cv_hu', realtype),
#		('cv_saved_tq5', realtype),
#		('cv_jcur', ctypes.c_int),
#		('cv_tolsf', realtype),
#		('cv_qmax_alloc', ctypes.c_int),
#		('cv_indx_acor', ctypes.c_int),
#		('cv_setupNonNull', ctypes.c_int),
#		('cv_VabstolMallocDone', ctypes.c_int),
#		('cv_MallocDone', ctypes.c_int),  
#		('cv_ehfun', ctypes.c_void_p),
#		('cv_eh_data', ctypes.c_void_p),
#		('cv_errfp', ctypes.c_void_p),
#		('cv_sldeton', ctypes.c_int),
#		('cv_ssdat', (realtype*6)*4),
#		('cv_nscon', ctypes.c_int),
#		('cv_nor', ctypes.c_long),
#		('cv_gfun', ctypes.c_void_p),
#		('cv_nrtfn', ctypes.c_int),
#		('cv_g_data', ctypes.c_void_p),
#		('*cv_iroots', ctypes.c_int),
#		('cv_tlo', realtype),
#		('cv_thi', realtype),
#		('cv_trout', realtype),
#		('*cv_glo', realtype),
#		('*cv_ghi', realtype),
#		('*cv_grout', realtype),
#		('cv_toutc', realtype),
#		('cv_ttol', realtype),
#		('cv_taskc', ctypes.c_int),
#		('cv_irfnd', ctypes.c_int),
#		('cv_nge', ctypes.c_long),
#	]

class CVodeMemObj(object):
	"""The CVodeMemObj class exists to provide automated memory management of the underlying void *cvodemem integrator memory, and should never be instantiated directly"""
	def __init__(self, obj):
		self.obj = obj
		self.dealloc = False

	def __del__(self):
		if self.dealloc:
			p = ctypes.c_void_p()
			p.value = self.obj
			cvode.CVodeFree(ctypes.byref(p))
cvode.CVodeFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvode.CVodeFree.restype = None

CVRhsFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVRhsFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a RHS function.

The callable python object must take exactly 4 parameters, which will be passed in as
	time_step (float)	the current value of the independent variable
	y (NVector)		the vector of current dependent values
	ydot (NVector)		undefined values, contents should be set to the new values of y
	f_data (c_void_p)	pointer to user data set by CVodeSetFdata

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVRhsFn)
	exec 'def __CallbackInterface_%s(t, y, ydot, f_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(ydot), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVRhsFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVRootFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype), ctypes.c_void_p)
def WrapCallbackCVRootFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a root finding function.

The callable python object must take exactly 4 parameters, which will be passed in as
	time_step (float)	the current value of the independent variable
	y (NVector)		the vector of current dependent values
	gout (NVector)		undefined values, contents should be set to the roots desired. If any element of gout is zero upon return, CVode will return early specifying the numbers of roots found, i.e. the number of 0 elements of gout
	g_data (c_void_p)	pointer to user data set by CVodeRootInit

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVRootFn)
	exec 'def __CallbackInterface_%s(t, y, gout, g_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), gout, g_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVRootFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVEwtFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVEwtFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify an error weight function.

The callable python object must take exactly 3 parameters, which will be passed in as
	y (NVector)		the vector of current dependent variables
	ewt (NVector)		undefined values, contents should be set to the error weights
	e_data (c_void_p)	pointer to user data set by CVodeSetEwtFn

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVEwtFn)
	exec 'def __CallbackInterface_%s(y, ewt, e_data):\n\treturn __ActualCallback[%i](nvecserial.NVector(y), nvecserial.NVector(ewt), e_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVEwtFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVErrHandlerFn = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_void_p)
def WrapCallbackCVErrHandlerFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify an error handler function.

The callable python object must take exactly 4 parameters, which will be passed in as
	error_code (int)	the code of the error
	module (string)		the name of the module producing the error
	function_name (string)	the function name in which the error appeared
	message (string)	a message describing the nature of the error
	eh_data (c_void_p)	a pointer to user data set by CVodeSetErrHandlerFn

and must have no return value."""
	if func == None:
		return ctypes.cast(None, CVErrHandlerFn)
	exec 'def __CallbackInterface_%s(error_code, module, function, msg, eh_data):\n\treturn __ActualCallback[%i](error_code, module, function, msg, eh_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVErrHandlerFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVodeCreate(lmm, iter):
	"""Sets up a CVODE internal integrator structure, and returns a handle it in the form of a CVodeMemObj
	lmm (int)	is the type of linear multistep method to be used. The legal values are CV_ADAMS and CV_BDF.
		CV_ADAMS	(Adams-Moulton) is recommended for non-stiff problems.
		CV_BDF		(Backward Differentiation Formula) is recommended for stiff problems.
	iter (int)	specifies whether functional or newton iteration will be used. Legal values are: 
		CV_FUNCTIONAL	does not require linear algebra 
		CV_NEWTON	requires the solution of linear systems. Requires the specification of a CVODE linear solver. (Recommended for stiff problems)."""
	obj = cvode.CVodeCreate(lmm, iter)
	if obj == None:
		raise AssertionError("SUNDIALS ERROR: CVodeCreate() failed - returned NULL pointer")
	return CVodeMemObj(obj)
cvode.CVodeCreate.argtypes = [ctypes.c_int, ctypes.c_int]
cvode.CVodeCreate.restype = ctypes.c_void_p

def CVodeSetErrHandlerFn(cvodememobj, func,  eh_data):
	"""Sets a user provided error handling function.
	func	is a python callable taking the parameters
		error_code (int)	the code of the error
		module (string)		the name of the module producing the error
		function_name (string)	the function name in which the error appeared
		message (string)	a message describing the nature of the error
		eh_data (c_void_p)	a pointer to user data as set by eh_data
	
	eh_data (c_void_p)	a pointer to any user data you would like passed through to your custom handler. The specified function will be passed eh_data as its eh_data parameter."""
	ret = cvode.CVodeSetErrHandlerFn(cvodememobj.obj, WrapCallbackCVErrHandlerFn(func), eh_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrHanderFn() failed with flag %i"%(ret))
cvode.CVodeSetErrHandlerFn.argtypes = [ctypes.c_void_p, CVErrHandlerFn, ctypes.c_void_p]
cvode.CVodeSetErrHandlerFn.restype = ctypes.c_int

def CVodeSetErrFile(cvodememobj, fileobj):
	"""Sets the file where all error messages and warnings will be written if the default error handling function is used.
	fileobj (file object)	must be open for writing"""
	ret = cvode.CVodeSetErrFile(cvodememobj.obj, sundials_core.fdopen(fileobj.fileno, fileobj.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrFile() failed with flag %i"%(ret))
cvode.CVodeSetErrFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvode.CVodeSetErrFile.restype = ctypes.c_int

def CVodeSetFdata(cvodememobj, f_data):
	"""Sets the pointer to any user data. If called, whenever the user's right hand side function is evaluated, it will be passed f_data as its f_data parameter."""
	ret = cvode.CVodeSetFdata(cvodememobj.obj, f_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetFdata() failed with flag %i"%(ret))
cvode.CVodeSetFdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvode.CVodeSetFdata.restype = ctypes.c_int

def CVodeSetEwtFn(cvodememobj, func, e_data):
	"""Sets a user provided error weight function.

	func	is a python callable taking the following parameters which must return 0 in the case of success, otherwise non-zero
		y (NVector)		the vector of current dependent variables
		ewt (NVector)		undefined values, contents should be set to the error weights
		e_data (c_void_p)	a pointer to user data as set by e_data

	e_data (c_void_p)	is a pointer to any user data you would like passed through to your custom function. The specified function will be passed e_data as its e_data parameter."""
	ret = cvode.CVodeSetEwtFn(cvodememobj.obj, WrapCallbackCVEwtFn(func), e_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetEwtFn() failed with flag %i"%(ret))
cvode.CVodeSetEwtFn.argtypes = [ctypes.c_void_p, CVEwtFn, ctypes.c_void_p]
cvode.CVodeSetEwtFn.restype = ctypes.c_int

def CVodeSetMaxOrd(cvodememobj, maxord):
	"""Sets the maximum lmm order to be used by the solver. Defaults to 12 for CV_ADAMS, and 5 for CV_BDF
	maxord (int)"""
	ret = cvode.CVodeSetMaxOrd(cvodememobj.obj, maxord)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxOrd() failed with flag %i"%(ret))
cvode.CVodeSetMaxOrd.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxOrd.restype = ctypes.c_int

def CVodeSetMaxNumSteps(cvodememobj, maxsteps):
	"""Sets the maximum number of internal steps to be taken by the solver in its attempt to reach tout.
	maxsteps (int)"""
	ret = cvode.CVodeSetMaxNumSteps(cvodememobj.obj, maxsteps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNumSteps() failed with flag %i"%(ret))
cvode.CVodeSetMaxNumSteps.argtypes = [ctypes.c_void_p, ctypes.c_long] 
cvode.CVodeSetMaxNumSteps.restype = ctypes.c_int

def CVodeSetMaxHnilWarns(cvodememobj, mxhnil):
	"""Sets the maximum number of warning messages issued by the solver that t+h==t on the next internal step. A value of -1 means no such warnings are issued. Default is 10.
	mxhnil (int)"""
	ret = cvode.CVodeSetMaxHnilWarns(cvodememobj.obj, mxhnil)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxHnilWarns() failed with flag %i"%(ret))
cvode.CVodeSetMaxHnilWarns.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxHnilWarns.restype = ctypes.c_int

def CVodeSetStabLimDet(cvodememobj, stldet):
	"""Turns stability limit detection on or off
	stldet (int)	0 = OFF, 1 = ON

When CV_BDF is used and order is 3 or greater, CVsldet is called to detect stability limit. If limit is detected, the order is reduced. Default is off"""
	ret = cvode.CVodeSetStabLimDet(cvodememobj.obj, stldet)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStabLimDet() failed with flag %i"%(ret))
cvode.CVodeSetStabLimDet.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetStabLimDet.restype = ctypes.c_int

def CVodeSetInitStep(cvodememobj, hin):
	"""Sets initial step size. CVODE estimates this value by default.
	hin (float)"""
	ret = cvode.CVodeSetInitStep(cvodememobj.obj, hin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetInitStep() failed with flag %i"%(ret))
cvode.CVodeSetInitStep.argtypes = [ctypes.c_void_p, realtype] 
cvode.CVodeSetInitStep.restype = ctypes.c_int

def CVodeSetMinStep(cvodememobj, hmin):
	"""Sets the absolute minimum of step size allowed. Default is 0.0
	hmin (float)"""
	ret = cvode.CVodeSetMinStep(cvodememobj.obj, hmin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMinStep() failed with flag %i"%(ret))
cvode.CVodeSetMinStep.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetMinStep.restype = ctypes.c_int

def CVodeSetMaxStep(cvodememobj, hmax):
	"""Sets the maximum absolute value of step size allowed. Default is infinity.
	hmax (float)"""
	ret = cvode.CVodeSetMaxStep(cvodememobj.obj, hmax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxStep() failed with flag %i"%(ret))
cvode.CVodeSetMaxStep.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetMaxStep.restype = ctypes.c_int

def CVodeSetStopTime(cvodememobj, tstop):
	"""Sets the independent variable value past which the solution is not to proceed. Default is infinity.
	tstop (float)"""
	ret = cvode.CVodeSetStopTime(cvodememobj.obj, tstop)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStopTime() failed with flag %i"%(ret))
cvode.CVodeSetStopTime.argtypes = [ctypes.c_void_p, realtype] 
cvode.CVodeSetStopTime.restype = ctypes.c_int

def CVodeSetMaxErrTestFails(cvodememobj, maxnef):
	"""Sets the maximum number of error test failures in attempting one step. Default is 7.
	maxnef (int)"""
	ret = cvode.CVodeSetMaxErrTestFails(cvodememobj.obj, maxnef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxErrTestFails() failed with flag %i"%(ret))
cvode.CVodeSetMaxErrTestFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxErrTestFails.restype = ctypes.c_int

def CVodeSetMaxNonlinIters(cvodememobj, maxcor):
	"""Sets the maximum number of nonlinear solver iterations at one solution. Default is 3.
	maxcor (int)"""
	ret = cvode.CVodeSetMaxNonlinIters(cvodememobj.obj, maxcor)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNonlinIters() failed with flag %i"%(ret))
cvode.CVodeSetMaxNonlinIters.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxNonlinIters.restype = ctypes.c_int

def CVodeSetMaxConvFails(cvodememobj, maxncf):
	"""Sets the maximum number of convergence failures allowed in attempting one step. Default is 10.
	maxncf (int)"""
	ret = cvode.CVodeSetMaxConvFails(cvodememobj.obj, maxncf)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxConvFails() failed with flag %i"%(ret))
cvode.CVodeSetMaxConvFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxConvFails.restype = ctypes.c_int

def CVodeSetNonlinConvCoef(cvodememobj, nlscoef):
	"""Sets the coefficient in the nonlinear convergence test. Default is 0.1
	nlscoef (float)"""
	ret = cvode.CVodeSetNonlinConvCoef(cvodememobj.obj, nlscoef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetNonlinConvCoef() failed with flag %i"%(ret))
cvode.CVodeSetNonlinConvCoef.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetNonlinConvCoef.restype = ctypes.c_int

def CVodeSetIterType(cvodememobj, iter):
	"""Changes the current nonlinear iteration type. Initially set by CVodeCreate().
	iter (int)	specifies whether functional or newton iteration will be used. Legal values are: 
		CV_FUNCTIONAL	does not require linear algebra 
		CV_NEWTON	requires the solution of linear systems. Requires the specification of a CVODE linear solver. (Recommended for stiff problems)."""
	ret = cvode.CVodeSetIterType(cvodememobj.obj, iter)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetIterType() failed with flag %i"%(ret))
cvode.CVodeSetIterType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetIterType.restype = ctypes.c_int

def CVodeSetTolerances(cvodememobj, itol, reltol, abstol):
	"""Changes the integration tolerances between calls to CVode(). Initially (re)set by CVodeMalloc()/CVodeReInit().
	itol (int)		the type of tolerances to be used. Legal values are:
		CV_SS	scalar relative and absolute tolerances
		CV_SV	scalar relative tolerance and a vector of absolute tolerances.
		CV_WF	indicates that the user will provide a function to evaluate the error weights. In this case, reltol and abstol are ignored, and should be None.
	reltol (float)		the relative tolerance scalar
	abstol (float/NVector)	absolute tolerance(s)

The parameters itol, reltol, and abstol define a vector of error weights, ewt, with components
	ewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = CV_SS), or
	ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = CV_SV).
This vector is used in all error and convergence tests, which use a weighted RMS norm on all error-like vectors v:
	WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
where N is the problem dimension."""
	if itol == CV_SS:
		if type(abstol) == realtype:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvode.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS, CV_SV or CV_WF")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetTolerances() failed with flag %i"%(ret))
cvode.CVodeSetTolerances.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype, ctypes.c_void_p]
cvode.CVodeSetTolerances.restype = ctypes.c_int

def CVodeMalloc(cvodememobj, func, t0, y0, itol, reltol, abstol):
	"""CVodeMalloc allocates and initializes memory for a problem to be solved by CVODE.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	func (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		time_step (float)	the current value of the independent variable
		y (NVector)		the vector of current dependent values
		ydot (NVector)		undefined values, contents should be set to the new values of y
		f_data (c_void_p)	pointer to user data set by CVodeSetFdata
	t0 (float)			the initial value of the independent variable.
	y0 (NVector)			the initial condition vector y(t0).
	itol (int)			the type of tolerances to be used. Legal values are:
		CV_SS	scalar relative and absolute tolerances
		CV_SV	scalar relative tolerance and a vector of absolute tolerances.
		CV_WF	indicates that the user will provide a function to evaluate the error weights. In this case, reltol and abstol are ignored, and should be None.
	reltol (float)			the relative tolerance scalar
	abstol (float/NVector)		absolute tolerance(s)

The parameters itol, reltol, and abstol define a vector of error weights, ewt, with components
	ewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = CV_SS), or
	ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = CV_SV).
This vector is used in all error and convergence tests, which use a weighted RMS norm on all error-like vectors v:
	WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
where N is the problem dimension."""
	if itol == CV_SS:
		if type(abstol) == realtype:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvode.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeMalloc() failed with flag %i"%(ret))
	cvodememobj.dealloc = True
cvode.CVodeMalloc.argtypes = [ctypes.c_void_p, CVRhsFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvode.CVodeMalloc.restype = ctypes.c_int

def CVodeReInit(cvodememobj, func, t0, y0, itol, reltol, abstol):
	"""CVodeReInit re-initializes CVode for the solution of a problem, where a prior call to CVodeMalloc has been made with the same problem size N. CVodeReInit performs the same input checking and initializations that CVodeMalloc does.  But it does no memory allocation, assuming that the existing internal memory is sufficient for the new problem.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	func (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		time_step (float)	the current value of the independent variable
		y (NVector)		the vector of current dependent values
		ydot (NVector)		undefined values, contents should be set to the new values of y
		f_data (c_void_p)	pointer to user data set by CVodeSetFdata
	t0 (float)			the initial value of the independent variable.
	y0 (NVector)			the initial condition vector y(t0).
	itol (int)			the type of tolerances to be used. Legal values are:
		CV_SS	scalar relative and absolute tolerances
		CV_SV	scalar relative tolerance and a vector of absolute tolerances.
		CV_WF	indicates that the user will provide a function to evaluate the error weights. In this case, reltol and abstol are ignored, and should be None.
	reltol (float)			the relative tolerance scalar
	abstol (float/NVector)		absolute tolerance(s)

The parameters itol, reltol, and abstol define a vector of error weights, ewt, with components
	ewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = CV_SS), or
	ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = CV_SV).
This vector is used in all error and convergence tests, which use a weighted RMS norm on all error-like vectors v:
	WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
where N is the problem dimension.

The use of CVodeReInit requires that the maximum method order, maxord, is no larger for the new problem than for the problem specified in the last call to CVodeMalloc. This condition is automatically fulfilled if the multistep method parameter lmm is unchanged (or changed from CV_ADAMS to CV_BDF) and the default value for maxord is specified."""
	if itol == CV_SS:
		if type(abstol) == realtype:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvode.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeReInit() failed with flag %i"%(ret))
cvode.CVodeReInit.argtypes = [ctypes.c_void_p, CVRhsFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvode.CVodeReInit.restype = ctypes.c_int

def CVodeRootInit(cvodememobj, nrtfn, func, g_data):
	"""CVodeRootInit initializes a rootfinding problem to be solved during the integration of the ODE system. It must be called after CVodeCreate, and before CVode. The arguments are:
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	nrtfn (int)			number of functions whose roots are to be found. nrtfn must be positive.
	func (function)			a python callable which will be automatically wrapped by WrapCallbackCVRootFn() which defines the functions g_i (0 < i < nrtfn) whose roots are sought. The function should take the the following parameters
		time_step (float)	the current value of the independent variable
		y (NVector)		the vector of current dependent values
		gout (NVector)		undefined values, contents should be set to the roots desired. If any element of gout is zero upon return, CVode will return early specifying the numbers of roots found, i.e. the number of 0 elements of gout
		g_data (c_void_p)	pointer to user data set by CVodeRootInit
	g_data (c_void_p)		a pointer to user data that will be passed to the user's root evaluation function every time it is called. The function receives the value of g_data in its g_data parameter

If a new problem is to be solved with a call to CVodeReInit, where the new problem has no root functions but the prior one did, then call CVodeRootInit with nrtfn = 0."""
	ret = cvode.CVodeRootInit(cvodememobj.obj, nrtfn, WrapCallbackCVRootFn(func), g_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeRootInit() failed with flag %i"%(ret))
cvode.CVodeRootInit.argtypes = [ctypes.c_void_p, ctypes.c_int, CVRootFn, ctypes.c_void_p]
cvode.CVodeRootInit.restype = ctypes.c_int

def CVode(cvodememobj, tout, yout, tret, itask):
	"""CVode integrates an ODE system over an interval in t. 
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	tout (float)			the next time at which a computed solution is desired.
	yout (NVector)			the computed solution vector. In CV_NORMAL mode with no errors and no roots found, yout = y(tout).
	tret (*realtype)		is set to the time reached by the solver
	itask (int)			is one of (See CVODE documentation for more details.)
		CV_NORMAL		the solver integrates from its current internal t value to a point at or beyond tout, then interpolates to t = tout and returns y(tout) in the user- allocated vector yout.
		CV_ONE_STEP		the solver takes one internal time step and returns in yout the value of y at the new internal time. In this case, tout is used only during the first call to CVode to determine the direction of integration and the rough scale of the t variable.
		CV_NORMAL_TSTOP		the solver returns the solution at tstop if that comes sooner than tout. The time reached by the solver is placed in tret. The user is responsible for allocating the memory for tret.
		CV_ONE_STEP_TSTOP	the solver returns the solution at tstop if that comes sooner than nect internal timestep. The time reached by the solver is placed in tret. The user is responsible for allocating the memory for tret.

CVode can return any of the following values:
	CV_SUCCESS:       CVode succeeded and no roots were found.  
	CV_ROOT_RETURN:   CVode succeeded, and found one or more roots. If nrtfn > 1, call CVodeGetRootInfo to see which g_i were found to have a root at (*tret).  
	CV_TSTOP_RETURN:  CVode succeeded and returned at tstop.  
	CV_MEM_NULL:      The cvode_mem argument was NULL.  
	CV_NO_MALLOC:     cvode_mem was not allocated.  
	CV_ILL_INPUT:     One of the inputs to CVode is illegal. This includes the situation when a component of the error weight vectors becomes < 0 during internal time-stepping.  It also includes the situation where a root of one of the root functions was found both at t0 and very near t0.  The ILL_INPUT flag will also be returned if the linear solver routine CV--- (called by the user after calling CVodeCreate) failed to set one of the linear solver-related fields in cvode_mem or if the linear solver's init routine failed. In any case, the user should see the printed error message for more details.  
	CV_TOO_MUCH_WORK: The solver took mxstep internal steps but could not reach tout. The default value for mxstep is MXSTEP_DEFAULT = 500.  
	CV_TOO_MUCH_ACC:  The solver could not satisfy the accuracy demanded by the user for some internal step.  
	CV_ERR_FAILURE:   Error test failures occurred too many times (= MXNEF = 7) during one internal time step or occurred with |h| = hmin.  
	CV_CONV_FAILURE:  Convergence test failures occurred too many times (= MXNCF = 10) during one internal time step or occurred with |h| = hmin.  
	CV_LINIT_FAIL:    The linear solver's initialization function failed.  
	CV_LSETUP_FAIL:   The linear solver's setup routine failed in an unrecoverable manner.  
	CV_LSOLVE_FAIL:   The linear solver's solve routine failed in an unrecoverable manner."""
	ret = cvode.CVode(cvodememobj.obj, tout, yout.data, tret, itask)
	return ret
cvode.CVode.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype), ctypes.c_int]
cvode.CVode.restype = ctypes.c_int

def CVodeGetDky(cvodememobj, t, k, dky):
	"""CVodeGetDky computes the kth derivative of the y function at time t.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	t (float)			the time at which the kth derivative of y is evaluated. The legal range for t is [tn-hu,tn] as described in the CVODE documentation, where
		tn	denotes the current internal time reached
		hu	denotes the last internal step size successfully used by the solver
	k (int)				the order of the derivative of y to be computed. The legal range for k is [0,qu] as described in the CVODE documentation.
		qu	denotes the order last used
	dky (NVector)			the kth derivative of the right hand side function will be will be placed in this vector; dky = [((d/dy)^k)y](t)

It is only legal to call this function after a successful return from CVode."""
	ret = cvode.CVodeGetDky(cvodememobj.obj, t, k, dky.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetDky() failed with flag %i"%(ret))
cvode.CVodeGetDky.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)] 
cvode.CVodeGetDky.restype = ctypes.c_int

def CVodeGetWorkSpace(cvodememobj):
	"""CVodeGetWorkSpace returns the CVODE real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrw = ctypes.c_long()
	leniw = ctypes.c_long()
	ret = cvode.CVodeGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrw), ctypes.byref(leniw))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetWorkSpace() failed with flag %i"%(ret))
	return (lenrw, leniw)
cvode.CVodeGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetWorkSpace.restype = ctypes.c_int

def CVodeGetNumSteps(cvodememobj):
	"""CVodeGetNumSteps returns the cumulative number of internal steps taken by the solver
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumSteps(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSteps() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumSteps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumSteps.restype = ctypes.c_int

def CVodeGetNumRhsEvals(cvodememobj):
	"""CVodeGetNumRhsEvals returns the number of calls to the user's f function
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumRhsEvals.restype = ctypes.c_int

def CVodeGetNumLinSolvSetups(cvodememobj):
	"""CVodeGetNumLinSolvSetups returns the number of calls made to the linear solver's setup routine
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumLinSolvSetups(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumLinSolvSetups() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumLinSolvSetups.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumLinSolvSetups.restype = ctypes.c_int

def CVodeGetNumErrTestFails(cvodememobj):
	"""CVodeGetNumErrTestFails returns the number of local error test failures that have occured
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumErrTestFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumErrTestFails() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumErrTestFails.restype = ctypes.c_int

def CVodeGetLastOrder(cvodememobj):
	"""CVodeGetLastOrder returns the order used during the last internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVodeGetLastOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastOrder() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetLastOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetLastOrder.restype = ctypes.c_int

def CVodeGetCurrentOrder(cvodememobj):
	"""CVodeGetCurrentOrder returns the order to be used on the next internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVodeGetCurrentOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentOrder() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetCurrentOrder.restype = ctypes.c_int

def CVodeGetNumStabLimOrderReds(cvodememobj):
	"""CVodeGetNumStabLimOrderReds returns the number of order reductions due to stability limit detection
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumStabLimOrderReds(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumStabLimOrderReds() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumStabLimOrderReds.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumStabLimOrderReds.restype = ctypes.c_int

def CVodeGetActualInitStep(cvodememobj):
	"""CVodeGetActualInitStep returns the actual initial step size used by CVODE
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetActualInitStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetActualInitStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetActualInitStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetActualInitStep.restype = ctypes.c_int

def CVodeGetLastStep(cvodememobj):
	"""CVodeGetLastStep returns the step size for the last internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetLastStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetLastStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetLastStep.restype = ctypes.c_int

def CVodeGetCurrentStep(cvodememobj):
	"""CVodeGetCurrentStep returns the step size to be attempted on the next internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetCurrentStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetCurrentStep.restype = ctypes.c_int

def CVodeGetCurrentTime(cvodememobj):
	"""CVodeGetCurrentTime returns the current internal time reached by the solver
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetCurrentTime(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentTime() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentTime.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetCurrentTime.restype = ctypes.c_int

def CVodeGetTolScaleFactor(cvodememobj):
	"""CVodeGetTolScaleFactor returns a suggested factor by which the user's tolerances should be scaled when too much accuracy has been requested for some internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetTolScaleFactor(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetTolScaleFactor() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetTolScaleFactor.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetTolScaleFactor.restype = ctypes.c_int

def CVodeGetErrWeights(cvodememobj, eweight):
	"""CVodeGetErrWeights returns the current error weight vector in eweight.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	eweight	(NVector)		is set to the current error weight vector"""
	if type(cvodememobj).__name__ == 'CVodeMemObj':
		cvo = cvodememobj.obj
	else:
		cvo = cvodememobj
	ret = cvode.CVodeGetErrWeights(cvo, eweight.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetErrWeights() failed with flag %i"%(ret))
cvode.CVodeGetErrWeights.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvode.CVodeGetErrWeights.restype = ctypes.c_int

def CVodeGetEstLocalErrors(cvodememobj, ele):
	"""CVodeGetEstLocalErrors returns the vector of estimated local errors in ele.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	ele (NVector)			is set to the vector of current estimated local errors"""
	ret = cvode.CVodeGetEstLocalErrors(cvodememobj.obj, ele.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetEstLocalErrors() failed with flag %i"%(ret))
cvode.CVodeGetEstLocalErrors.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvode.CVodeGetEstLocalErrors.restype = ctypes.c_int

def CVodeGetNumGEvals(cvodememobj):
	"""CVodeGetNumGEvals returns the number of calls to the user's g function (for rootfinding)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumGEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumGEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumGEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumGEvals.restype = ctypes.c_int

def CVodeGetRootInfo(cvodememobj, numroots):
	"""CVodeGetRootInfo returns the indices for which the function g_i was found to have a root. A list of zeroes and ones is returned, where a one in index postion i indicates a root for g_i was found.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	numroots (int)			the number of functions for which roots should be found."""
	rootsfound = (ctypes.c_int*numroots)()
	ret = cvode.CVodeGetRootInfo(cvodememobj.obj, rootsfound)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVode() failed with flag %i"%(ret))
	return [rootsfound[i] for i in range(numroots)]
cvode.CVodeGetRootInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetRootInfo.restype = ctypes.c_int
	
def CVodeGetIntegratorStats(cvodememobj):
	"""A convenience function which returns a tuple of all available integrator stats.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()

Returns a tuple of the format:
 
	(nsteps, nfevals, nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur) where: 
		nsteps	 	cumulative number of internal steps taken by the solver 
		nfevals	 	number of rhs evaluations 
		nlinsetups	number of calls to the linear setup function 
		netfails	number of local error test failures 
		qlast	 	order used during ast internal step 
		qcur	 	order to be used on next internal step 
		hinused	 	actual initial step size ised by CVODE 
		hlast	 	last step size used by CVODE 
		hcur	 	next step size to be used by CVODE 
		tcur		current time"""
	nsteps = ctypes.c_long()
	nfevals = ctypes.c_long()
	nlinsetups = ctypes.c_long()
	netfails = ctypes.c_long()
	qlast = ctypes.c_int()
	qcur = ctypes.c_int()
	hinused = realtype()
	hlast = realtype()
	hcur = realtype()
	tcur = realtype()
	ret = cvode.CVodeGetIntegratorStats(cvodememobj.obj, ctypes.byref(nsteps), ctypes.byref(nfevals), ctypes.byref(nlinsetups), ctypes.byref(netfails), ctypes.byref(qlast), ctypes.byref(qcur), ctypes.byref(hinused), ctypes.byref(hlast), ctypes.byref(hcur), ctypes.byref(tcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetIntegratorStats() failed with flag %i"%(ret))
	return (nsteps.value, nfevals.value, nlinsetups.value, netfails.value, qlast.value, qcur.value, hinused.value, hlast.value, hcur.value, tcur.value)
cvode.CVodeGetIntegratorStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(realtype), ctypes.POINTER(realtype), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
cvode.CVodeGetIntegratorStats.restype = ctypes.c_int

def CVodeGetNumNonlinSolvIters(cvodememobj):
	"""CVodeGetNumNonlinSolvIters returns the number of nonlinear solver iterations performed.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumNonlinSolvIters(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvIters() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumNonlinSolvIters.restype = ctypes.c_int

def CVodeGetNumNonlinSolvConvFails(cvodememobj):
	"""CVodeGetNumNonlinSolvConvFails returns the number of nonlinear convergence failures.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumNonlinSolvConvFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvConvFails() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumNonlinSolvConvFails.restype = ctypes.c_int

def CVodeGetNonlinSolvStats(cvodememobj):
	"""A convenience function that provides the nonlinear solver optional outputs in a tuple (NumNonlinSolvIters, NumNonlinsolvConvFails).
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nniters = ctypes.c_long()
	nncfails = ctypes.c_long()
	ret = cvode.CVodeGetNonlinSolvStats(cvodememobj.obj, ctypes.byref(nniters), ctypes.byref(nncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNonlinSolvStats() failed with flag %i"%(ret))
	return (nniters.value, nncfails.value)
cvode.CVodeGetNonlinSolvStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNonlinSolvStats.restype = ctypes.c_int

def CVodeGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVODE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVodeGetReturnFlagName(flag)
cvode.CVodeGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVodeGetReturnFlagName.restype = ctypes.c_char_p

########################
# sundials_iterative.h #
########################

PREC_NONE = 0
PREC_LEFT = 1
PREC_RIGHT = 2
PREC_BOTH = 3

MODIFIED_GS = 1
CLASSICAL_GS = 2

ATimesFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackATimesFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions requiring an ATimesFn function.

The callable python object must take exactly 4 parameters, which will be passed in as
	A_data (c_voip)	a pointer to user data as passed in by SpbcgSolve, SpgmrSolve, or SptfqmrSolve
	v (NVector)	the vector to be multiplied by A
	z (NVector)	undefined values, content should modified to store the result of A times v

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if (func == None):
		return ctypes.cast(None, ATimesFn)
	exec 'def __CallbackInterface_%s(A_data, v, z):\n\treturn __ActualCallback[%i](A_data, nvecserial.NVector(v), nvecserial.NVector(z))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = ATimesFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

PSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int)
def WrapCallbackPSolveFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions requiring an PSolve function.

A PSolveFn solves the preconditioner equation Pz = r for the vector z.
The callable python object must take exactly 4 parameters, which will be passed in as
	P_data (c_voip)	a pointer to user data as passed in by SpbcgSolve, SpgmrSolve, or SptfqmrSolve
	r (NVector)	the vector for which the preconditioner equation Pz = r must be solves
	z (NVector)	undefined values, content should modified to store the solution
	lr (int)	1 for left preconditioning, 2 for right. If preconditioning is on one side only, lr can be ignored.

and must return an integer of 0 in the case of no error. Otherwise a negative value indicates an unrecoverable error, and a positive value, a recoverable one, in which the calling routine may reattempt the solution after updating preconditioner data."""
	if (func == None):
		return ctypes.cast(None, PSolveFn)
	exec 'def __CallbackInterface_%s(P_data, r, z, lr):\n\treturn __ActualCallback[%i](P_data, nvecserial.NVector(r), nvecserial.NVector(z), lr)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = PSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def ModifiedGS(v, h, k, p, new_vk_norm):
 	"""ModifiedGS performs a modified Gram-Schmidt orthogonalization of the N_Vector v[k] against the p unit N_Vectors at v[k-1], v[k-2], ..., v[k-p].

	v (NVectorArray)	an array of (k+1) N_Vectors v[i], i=0, 1, ..., k.  v[k-1], v[k-2], ..., v[k-p] are assumed to have L2-norm equal to 1.
	h (**realtype)		the output k by k Hessenberg matrix of inner products. This matrix must be allocated row-wise so that the (i,j)th entry is h[i][j]. The inner products (v[i],v[k]), i=i0, i0+1, ..., k-1, are stored at h[i][k-1]. Here i0=MAX(0,k-p).
	k (int)			the index of the vector in the v array that needs to be orthogonalized against previous vectors in the v array.
	p (int)			the number of previous vectors in the v array against which v[k] is to be orthogonalized.
	new_vk_norm (*realtype)	a pointer to memory allocated by the caller to hold the Euclidean norm of the orthogonalized vector v[k].

If (k-p) < 0, then ModifiedGS uses p=k. The orthogonalized v[k] is NOT normalized and is stored over the old v[k]. Once the orthogonalization has been performed, the Euclidean norm of v[k] is stored in (*new_vk_norm).

ModifiedGS returns 0 to indicate success. It cannot fail."""
	ret = cvode.ModifiedGS(v.data, h, k, p, new_vk_norm)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ModifiedGS() failed with flag %i"%(ret))
cvode.ModifiedGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype)]
cvode.ModifiedGS.restype = ctypes.c_int

def ClassicalGS(v, h, k, p, new_vk_norm, temp, s):
	"""ClassicalGS performs a classical Gram-Schmidt orthogonalization of the N_Vector v[k] against the p unit N_Vectors at v[k-1], v[k-2], ..., v[k-p].

	v (NVectorArray)	an array of (k+1) N_Vectors v[i], i=0, 1, ..., k.  v[k-1], v[k-2], ..., v[k-p] are assumed to have L2-norm equal to 1.
	h (**realtype)		the output k by k Hessenberg matrix of inner products. This matrix must be allocated row-wise so that the (i,j)th entry is h[i][j]. The inner products (v[i],v[k]), i=i0, i0+1, ..., k-1, are stored at h[i][k-1]. Here i0=MAX(0,k-p).
	k (int)			the index of the vector in the v array that needs to be orthogonalized against previous vectors in the v array.
	p (int)			the number of previous vectors in the v array against which v[k] is to be orthogonalized.
	new_vk_norm (*realtype)	a pointer to memory allocated by the caller to hold the Euclidean norm of the orthogonalized vector v[k].
	temp (NVector)		which can be used as workspace by the ClassicalGS routine.
	s (*realtype)		a length k array of realtype which can be used as workspace by the ClassicalGS routine.

ClassicalGS returns 0 to indicate success. It cannot fail."""
	ret = cvode.ClassicalGS(v.data, h, k, p, new_vk_norm, temp.data, s)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ClassicalGS() failed with flag %i"%(ret))
cvode.ClassicalGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype)]
cvode.ClassicalGS.restype = ctypes.c_int

def QRfact(n, h, q, job):
	"""QRfact performs a QR factorization of the Hessenberg matrix H.
	n (int)		the problem size; the matrix H is (n+1) by n.
	h (**realtype)	the (n+1) by n Hessenberg matrix H to be factored. It is stored row-wise.
	q (*realtype)	an array of length 2*n containing the Givens rotations computed by this function. A Givens rotation has the form:
				| c  -s |
				| s   c |.
			The components of the Givens rotations are stored in q as (c, s, c, s, ..., c, s).
	job (int)	a control flag. If job==0, then a new QR factorization is performed. If job!=0, then it is assumed that the first n-1 columns of h have already been factored and only the last column needs to be updated.

QRfact returns 0 if successful. If a zero is encountered on the diagonal of the triangular factor R, then QRfact returns the equation number of the zero entry, where the equations are numbered from 1, not 0. If QRsol is subsequently called in this situation, it will return an error because it could not divide by the zero diagonal entry."""
	ret = cvode.QRfact(n, h, q, job)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRfact() failed with flag %i"%(ret))
cvode.QRfact.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.c_int]
cvode.QRfact.restype = ctypes.c_int

def QRsol(n, h, q, b):
	"""QRsol solves the linear least squares problem

min (b - H*x, b - H*x), x in R^n,

where H is a Hessenberg matrix, and b is in R^(n+1).

It uses the QR factors of H computed by QRfact.
	n (int)		the problem size; the matrix H is (n+1) by n.
	h (**realtype)	a matrix (computed by QRfact) containing the upper triangular factor R of the original Hessenberg matrix H.
	q (*realtype)	an array of length 2*n (computed by QRfact) containing the Givens rotations used to factor H. A Givens rotation has the form:
				| c  -s |
				| s   c |.
			The components of the Givens rotations are stored in q as (c, s, c, s, ..., c, s).
	b (*realype)	the (n+1)-vector appearing in the least squares problem above.

On return, b contains the solution x of the least squares problem, if QRsol was successful.

QRsol returns a 0 if successful.  Otherwise, a zero was encountered on the diagonal of the triangular factor R.  In this case, QRsol returns the equation number (numbered from 1, not 0) of the zero entry."""
	ret = cvode.QRsol(n, h, q, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRsol() failed with flag %i"%(ret))
cvode.QRsol.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
cvode.QRsol.restype = ctypes.c_int


#########################
# sundials_smalldense.h #
#########################

def denalloc(m, n):
	"""Allocates storage for an m by n dense matrix and returns a pointer to the newly allocated storage if successful. If the memory request cannot be satisfied, a ValueError Exception is raised.

The underlying type of the dense matrix returned is a double pointer to realtype. If we allocate a dense matrix a by a = denalloc(m, n), then a[j][i] references the (i,j)-th element of the matrix a, 0 <= i < m, 0 <= j < n,  and a[j] is a pointer to the first element in the jth column of a. The location a[0] contains a pointer to m*n contiguous locations which contain the elements of a."""
	ret = cvode.denalloc(m, n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denalloc could not allocate memory')
cvode.denalloc.argtypes = [ctypes.c_long, ctypes.c_long]
cvode.denalloc.restype = ctypes.POINTER(ctypes.POINTER(realtype))

def denallocpiv(n):
	"""Allocates an array of n long int, and returns a pointer to the first element in the array if successful, otherwise a ValueError exception is raised."""
	ret = cvode.denallocpiv(n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denallocpiv could not allocate memory')
cvode.denallocpiv.argtypes = [ctypes.c_long]
cvode.denallocpiv.restype = ctypes.POINTER(ctypes.c_long)

def denGETRF(a, m, n, p):
	"""Factors the m by n dense matrix a, m>=n. It overwrites the elements of a with its LU factors and keeps track of the pivot rows chosen in the pivot array p. A successful LU factorization leaves the matrix a and the pivot array p with the following information:
	(1) p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., n-1.
	(2) If the unique LU factorization of a is given by Pa = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1.0 on the diagonal, and U is an upper triangular matrix, then the upper triangular part of a (including its diagonal) contains U and the strictly lower trapezoidal part of a contains the multipliers, I-L.

Note that for square matrices (m=n), L is unit lower triangular.

denGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	return cvode.denGETRF(a, m, n, p)
cvode.denGETRF.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]
cvode.denGETRF.restype = ctypes.c_long

def denGETRS(a, n, p, b):
	"""Solves the n by n linear system a*x = b. It assumes that a has been LU factored and the pivot array p has been set by a successful call to denGETRF(a,n,n,p).
	
denGETRS does not check whether a is square! The solution x is written into the b array."""
	return cvode.denGETRS(a, n, p, b)
cvode.denGETRS.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
cvode.denGETRS.restype = None

def denzero(a, m, n):
	"""Sets all the elements of the m by n dense matrix a to be 0.0."""
	cvode.denzero(a, m, n)
cvode.denzero.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.denzero.restype = None

def dencopy(a, b, m, n):
	"""Copies the m by n dense matrix a into the m by n dense matrix b."""
	cvode.dencopy(a, b, m, n)
cvode.dencopy.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.dencopy.restype = None

def denscale(c, a, m, n):
	"""Scales every element in the m by n dense matrix a by c."""
	cvode.denscale(c, a, m, n)
cvode.denscale.argtypes = [realtype, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.denscale.restype = None

def denaddI(a, n):
	"""Increments the diagonal elements of the dense m by n matrix a by 1.0. (a_ii <= a_ii + 1, i=1,2,..n-1.)
	
denaddI is typically used with square matrices.
denaddI does NOT check for m >= n! Therefore, a segmentation fault will occur if m<n!
"""
	cvode.denaddI(a, n)
cvode.denaddI.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long]
cvode.denaddI.restype = None

def denfreepiv(p):
	"""Frees the pivot array p allocated by denallocpiv."""
	cvode.denfreepiv(p)
cvode.denfreepiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvode.denfreepiv.restype = None

def denfree(a):
	"""Frees the dense matrix a allocated by denalloc."""
	cvode.denfree(a)
cvode.denfree.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype))]
cvode.denfree.restype = None

def denprint(a, m, n):
	"""Prints the m by n dense matrix a to standard output as it would normally appear on paper. It is intended as a debugging tool with small values of m and n. The elements are printed using the %g/%lg/%Lg option. A blank line is printed before and after the matrix."""
	cvode.denprint(a, m, n)
cvode.denprint.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.denprint.restype = None

#####################
# sundials_spbcgs.h #
#####################

SPBCG_SUCCESS = 0
SPBCG_RES_REDUCED = 1
SPBCG_CONV_FAIL = 2
SPBCG_PSOLVE_FAIL_REC = 3
SPBCG_ATIMES_FAIL_REC = 4
SPBCG_PSET_FAIL_REC = 5

SPBCG_MEM_NULL = -1
SPBCG_ATIMES_FAIL_UNREC = -2
SPBCG_PSOLVE_FAIL_UNREC = -3
SPBCG_PSET_FAIL_UNREC = -4

class SpbcgMemRec(ctypes.Structure):
	"""Contains numerous fields that must be accessed by the SPBCG linear solver module.
	l_max (int)		maximum Krylov subspace dimension that SpbcgSolve will be permitted to use
	r (NVector) 		holds the scaled, preconditioned linear system residual
	r_star (NVector) 	holds the initial scaled, preconditioned linear system residual
	p (NVector)		used for workspace by the SPBCG algorithm
	q (NVector)		used for workspace by the SPBCG algorithm
	u (NVector)		used for workspace by the SPBCG algorithm
	Ap (NVector)		used for workspace by the SPBCG algorithm
	vtemp (NVector)		scratch vector used as temporary vector storage"""
	_fields_ = [
		("l_max", ctypes.c_int),
		("r_star", ctypes.POINTER(nvecserial._NVector)),
		("r", ctypes.POINTER(nvecserial._NVector)),
		("p", ctypes.POINTER(nvecserial._NVector)),
		("q", ctypes.POINTER(nvecserial._NVector)),
		("u", ctypes.POINTER(nvecserial._NVector)),
		("Ap", ctypes.POINTER(nvecserial._NVector)),
		("vtemp", ctypes.POINTER(nvecserial._NVector))
	]

SpbcgMem = ctypes.POINTER(SpbcgMemRec)

def SpbcgMalloc(l_max, vec_tmpl):
	"""Allocates additional memory needed by the SPBCG linear solver module.
	l_max (int)		maximum Krylov subspace dimension that SpbcgSolve will be permitted to use
	vec_tmpl (NVector)	implementation-specific template vector

If successful, SpbcgMalloc returns an SpbcgMemRec. 
an error occurs, then a NULL pointer is returned."""
	return cvode.SpbcgMalloc(l_max, vec_tmpl.data)
cvode.SpbcgMalloc.argtypes = [ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvode.SpbcgMalloc.restype = SpbcgMem

def SpbcgSolve(mem, A_data, x, b, pretype, delta, P_data, sx, sb, atimes, psolve, res_norm, nli, nps):
	"""Solves the linear system Ax = b by means of a scaled preconditioned Bi-CGSTAB (SPBCG) iterative method.
	mem (SpbcgMem)		handle returned by a previous call SpbcgMalloc
	A_data (c_void_p)	pointer to a data structure containing information about the coefficient matrix A (passed to user-supplied function referenced by atimes (function pointer))
	x (NVector)		contains initial guess x_0 upon entry, but which upon return contains an approximate solution of the linear system Ax = b (solution only valid if return value is either SPBCG_SUCCESS or SPBCG_RES_REDUCED)
	b (NVector)		is set to the right-hand side vector b of the linear system (undisturbed by function)
	pretype (int)		indicates the type of preconditioning to be used (see sundials_iterative.h)
	delta (realtype)	tolerance on the L2 norm of the scaled, preconditioned residual (if return value == SPBCG_SUCCESS, then ||sb*P1_inv*(b-Ax)||_L2 <= delta)
	P_data (c_void_p)	pointer to a data structure containing preconditioner information (passed to user-supplied function referenced by psolve (function pointer))
	sx (NVector)		contains positive scaling factors for x (pass sx = None if scaling NOT required)
	sb (NVector)		contains positive scaling factors for b (pass sb = None if scaling NOT required)
	atimes (func)		user-supplied routine responsible for computing the matrix-vector product Ax (see sundials_iterative.h)
	psolve (func)		user-supplied routine responsible for solving the preconditioned linear system Pz = r (ignored if pretype == PREC_NONE) (see sundials_iterative.h)
	res_norm (*realtype)	the L2 norm of the scaled, preconditioned residual (if return value is either SPBCG_SUCCESS or SPBCG_RES_REDUCED, then *res_norm = ||sb*P1_inv*(b-Ax)||_L2, where x is the computed approximate solution, sb is the diagonal scaling matrix for the right-hand side b, and P1_inv is the inverse of the left-preconditioner matrix)
	nli (*int)		is set to the total number of linear iterations performed
	nps (*int)		is set to the total number of calls made to the psolve routine"""
	ret = cvode.SpbcgSolve(mem, A_data, x.data, b.data, pretype, delta, P_data, sx.data, sb.data, WrapCallbackATimesFn(atimes), WrapCallbackPSolveFn(psolve), res_norm, nli, nps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgSolve() failed with flag %i"%(ret))
cvode.SpbcgSolve.argtypes = [SpbcgMem, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ATimesFn, PSolveFn, ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
cvode.SpbcgSolve.restype = ctypes.c_int

def SpbcgFree(mem):
	"""Deallocates any memory allocated implicitly by Spbcg 
	ret = cvode.SpbcgFree(mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgFree() failed with flag %i"%(ret))
cvode.SpbcgFree.argtypes = [SpbcgMem]
cvode.SpbcgFree.restype = None


################
# cvode_diag.h #
################

#CVDIAG return values
CVDIAG_SUCCESS = 0
CVDIAG_MEM_NULL = -1
CVDIAG_LMEM_NULL = -2
CVDIAG_ILL_INPUT = -3
CVDIAG_MEM_FAIL = -4

#Addtional last_flag values
CVDIAG_INV_FAIL = -5
CVDIAG_RHSFUNC_UNRECVR = -6
CVDIAG_RHSFUNC_RECVR = -7

def CVDiag(cvodememobj):
	"""Links the main integrator with the CVDIAG linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	ret = cvode.CVDiag(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiag() failed with flag %i"%(ret))
cvode.CVDiag.argtypes = [ctypes.c_void_p]
cvode.CVDiag.restype = ctypes.c_int

def CVDiagGetWorkSpace(cvodememobj):
	"""Returns the CVDIAG real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVDiagGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVDiagGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVDiagGetWorkSpace.restype = ctypes.c_intw
def CVDiagGetNumRhsEvals(cvodememobj):
	"""Returns the number of calls to the user f routine due to finite difference Jacobian evaluation.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()

Note: The number of diagonal approximate Jacobians formed is equal to the number of CVDiagSetup calls. This number is available through CVodeGetNumLinSolvSetups."""
	retval = ctypes.c_long(0)
	ret = cvode.CVDiagGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDiagGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVDiagGetNumRhsEvals.restype = ctypes.c_int

def CVDiagGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVDIAG interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int(0)
	ret = cvode.CVDiagGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetLastFlag() failed with flag %i"%(ret))
	return retval
cvode.CVDiagGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVDiagGetLastFlag.restype = ctypes.c_int

def CVDiagGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDIAG return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVDiagGetReturnFlagName(flag)
cvode.CVDiagGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVDiagGetReturnFlagName.restype = ctypes.c_char_p

##################
# cvode_bbdpre.h #
##################

CVBBDPRE_SUCCESS = 0
CVBBDPRE_PDATA_NULL = -11
CVBBDPRE_FUNC_UNRECVR = -12

CVLocalFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVLocalFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBBDPrecAlloc and CVBBDPrecReInit.

The callable python object must take exactly 5 parameters, which will be passed in as
	NLocal (ing)		the local vector size
	t (realtype)		the value of the independent variable value
	y (NVector)		the local real dependent variable vector
	g (NVector)		undefined values, must be set to the computation g(t,y) where g is an approximation of the RHS function (the case of g being mathematially identical to the RHS is allowed)
	g_data (c_void_p)	a pointer to the user-defined data block g_data, as set by a previous call to CVodeSetFdata

and must return an integer of 0 in the case of no error, a positive integer if a recovorable error occured, and a negative value in the case of an unrecoverable error."""
	if func == None:
		return ctypes.cast(None, CVLocalFn)
	exec 'def __CallbackInterface_%s(Nlocal, t, y, g, f_data):\n\treturn __ActualCallback[%i](Nlocal, t, nvecserial.NVector(y), nvecserial.NVector(g), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVLocalFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVCommFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVCommFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBBDPrecAlloc and CVBBDPrecReInit.

The callable python object performs all interprocess communication necessary to evaluate the approximate right-hand side function as described in WrapCallbackCVLocalFn, and must take exactly 4 parameters, which will be passed in as
	NLocal (ing)		the local vector size
	t (realtype)		the value of the independent variable value
	y (NVector)		the local real dependent variable vector
	f_data (c_void_p)	a pointer to the user-defined data block f_data, as set by a previous call to CVodeSetFdata

and must return an integer of 0 in the case of no error, a positive integer if a recovorable error occured, and a negative value in the case of an unrecoverable error.

The CVCommFn cfn is expected to save communicated data in space defined within the structure f_data. Each call to the CVCommFn function is preceded by a call to the CVRhsFn f with the same (t,y) arguments. Thus cfn can omit any communications done by f if relevant to the evaluation of g. If all necessary communication was done by f, the user can pass None for cfn in CVBBDPrecAlloc/CVBBEDPrecReInit."""
	if func == None:
		return ctypes.cast(None, CVCommFn)
	exec 'def __CallbackInterface_%s(Nlocal, t, y, f_data):\n\treturn __ActualCallback[%i](Nlocal, t, nvecserial.NVector(y), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVCommFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVBBDPrecAlloc(cvodememobj, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn):
	"""CVBBDPrecAlloc allocates and initializes a CVBBDData structure to be passed to CVSp* (and used by CVBBDPrecSetup and CVBBDPrecSolve).
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate() 
	Nlocal(int)			the length of the local block of the vectors y etc. on the current processor.  
	mudq (int)			the upper half-bandwidths to be used in the difference quotient computation of the local Jacobian block.  
	mldq (int)			the lower half-bandwidths to be used in the difference quotient computation of the local Jacobian block.  
	mukeep (int)			the upper half-bandwidths of the retained banded approximation to the local Jacobian block.  
	mlkeep (int)			the lower half-bandwidths of the retained banded approximation to the local Jacobian block.  
	dqrely (realtype)		the relative increment in components of y used in the difference quotient approximations. To specify the default, pass 0.  The default is dqrely = sqrt(unit roundoff).
	gloc (func)			the name of the user-supplied function g(t,y) that approximates f and whose local Jacobian blocks are to form the preconditioner.  
	cfn (func)			the name of the user-defined function that performs necessary interprocess communication for the execution of gloc.  

CVBBDPrecAlloc returns an object reperesenting the storage allocated."""
	ret = cvode.CVBBDPrecAlloc(cvodememobj.obj, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, WrapCallbackCVLocalFn(gloc), WrapCallbackCVCommFn(cfn))
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecAlloc() failed to allocate memory")
	return ret
cvode.CVBBDPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, CVLocalFn, CVCommFn]
cvode.CVBBDPrecAlloc.restype = ctypes.c_void_p

def CVBBDSptfqmr(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSptfqmr links the CVBBDPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions: 
	1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory 
	2) Sets the preconditioner data structure for CVSPTFQMR 
	3) Sets the preconditioner setup routine for CVSPTFQMR 
	4) Sets the preconditioner solve routine for CVSPTFQMR 

Its first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSptfqmr.  

Possible return values are: 
	CVSPILS_SUCCESS      if successful 
	CVSPILS_MEM_NULL     if the cvode memory was NULL 
	CVSPILS_LMEM_NULL    if the cvsptfqmr memory was NULL 
	CVSPILS_MEM_FAIL     if there was a memory allocation failure 
	CVSPILS_ILL_INPUT    if a required vector operation is missing 
	CVBBDPRE_PDATA_NULL  if the bbd_data was NULL"""
	ret = cvode.CVBBDSptfqmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSptfqmr() failed with flag %i"%(ret))
cvode.CVBBDSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSptfqmr.restype = ctypes.c_int

def CVBBDSpbcg(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSptfqmr links the CVBBDPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions:
	1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory;
	2) Sets the preconditioner data structure for CVSPTFQMR
	3) Sets the preconditioner setup routine for CVSPTFQMR
	4) Sets the preconditioner solve routine for CVSPTFQMR

Its first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSptfqmr."""
	ret = cvode.CVBBDSpbcg(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpbcg() failed with flag %i"%(ret))
cvode.CVBBDSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSpbcg.restype = ctypes.c_int

def CVBBDSpgmr(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSpbcg links the CVBBDPRE preconditioner to the CVSPBCG linear solver. It performs the following actions:
	1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory;
	2) Sets the preconditioner data structure for CVSPBCG
	3) Sets the preconditioner setup routine for CVSPBCG
	4) Sets the preconditioner solve routine for CVSPBCG

Its first 3 arguments are the same as for CVSpbcg (see cvspbcg.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSpbcg."""
	ret = cvode.CVBBDSpgmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpgmr() failed with flag %i"%(ret))
cvode.CVBBDSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSpgmr.restype = ctypes.c_int

def CVBBDPrecReInit(bbd_data, mudq, mldq, dqrely, gloc, cfn):
	"""CVBBDPrecReInit re-initializes the BBDPRE module when solving a sequence of problems of the same size with CVSPGMR/CVBBDPRE or CVSPBCG/CVBBDPRE or CVSPTFQMR/CVBBDPRE provided there is no change in Nlocal, mukeep, or mlkeep. After solving one problem, and after calling CVodeReInit to re-initialize the integrator for a subsequent problem, call CVBBDPrecReInit. Then call CVSpgmrSet* or CVSpbcgSet* or CVSptfqmrSet* functions if necessary for any changes to CVSPGMR, CVSPBCG, or CVSPTFQMR parameters, before calling CVode.

The first argument to CVBBDPrecReInit must be the pointer pdata that was returned by CVBBDPrecAlloc. All other arguments have the same names and meanings as those of CVBBDPrecAlloc."""
	ret = cvode.CVBBDPrecReInit(bbd_data, mudq, mldq, dqrely, WrapCallbackCVLocalFn(gloc), WrapCallbackCVCommFn(cfn))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecReInit() failed with flag %i"%(ret))
cvode.CVBBDPrecReInit.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, realtype, CVLocalFn, CVCommFn]
cvode.CVBBDPrecReInit.restype = ctypes.c_int

def CVBBDPrecFree(bbd_data):
	"""CVBBDPrecFree frees the memory block bbd_data allocated by the call to CVBBDAlloc."""
	ret = cvode.CVBBDPrecFree(bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecFree() failed with flag %i"%(ret))
cvode.CVBBDPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvode.CVBBDPrecFree.restype = None

def CVBBDPrecGetWorkSpace(bbd_data):
	"""CVBBDPrecGetWorkSpace returns the CVBBDPrec real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVBBDPrecGetWorkSpace(bbd_data, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVBBDPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVBBDPrecGetWorkSpace.restype = ctypes.c_int

def CVBBDPrecGetNumGfnEvals(bbd_data):
 	"""CVBBDPrecGetNumGfnEvals returns the number of calls to gfn."""
	ngevalsBBDP = ctypes.c_long(0)
	ret = cvode.CVBBDPrecGetNumGfnEvals(bbd_data, ctypes.byref(ngevalsBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecGetNumGfnEvals() failed with flag %i"%(ret))
	return ngevalsBBDP.value
cvode.CVBBDPrecGetNumGfnEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVBBDPrecGetNumGfnEvals.restype = ctypes.c_int

def CVBBDPrecGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBBDPRE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVBBDPrecGetReturnFlagName(flag)
cvode.CVBBDPrecGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVBBDPrecGetReturnFlagName.restype = ctypes.c_char_p

###################
# cvode_bandpre.h #
###################

CVBANDPRE_SUCCESS = 0
CVBANDPRE_PDATA_NULL = -11
CVBANDPRE_RHSFUNC_UNRECVR = -12

def CVBandPrecAlloc(cvodememobj, N, mu, ml):
	"""CVBandPrecAlloc allocates and initializes a CVBandPrecData structure to be passed to CVSp* (and subsequently used by CVBandPrecSetup and CVBandPrecSolve). The parameters of CVBandPrecAlloc are as follows:
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate() 
	N	the problem size.  
	mu	the upper half bandwidth.  
	ml	the lower half bandwidth.
CVBandPrecAlloc returns the storage pointer of type CVBandPrecData, or NULL if the request for storage cannot be satisfied.  
NOTE:	The band preconditioner assumes a serial implementation of the NVECTOR package. Therefore, CVBandPrecAlloc will first test for a compatible N_Vector internal representation by checking for required functions."""
	ret = cvode.CVBandPrecAlloc(cvodememobj.obj, N, mu, ml)
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecAlloc() failed to allocate memory")
	return ret
cvode.CVBandPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvode.CVBandPrecAlloc.restype = ctypes.c_void_p

def CVBPSptfqmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSptfqmr links the CVBANDPPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions: 
	1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPTFQMR 
	3) Sets the preconditioner setup routine for CVSPTFQMR 
	4) Sets the preconditioner solve routine for CVSPTFQMR 
Its first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSptfqmr."""
	ret = cvode.CVBPSptfqmr(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSptfqmr() failed with flag %i"%(ret))
cvode.CVBPSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBPSptfqmr.restype = ctypes.c_int

def CVBPSpbcg(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpbcg links the CVBANDPPRE preconditioner to the CVSPBCG linear solver. It performs the following actions: 
	1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPBCG 
	3) Sets the preconditioner setup routine for CVSPBCG 
	4) Sets the preconditioner solve routine for CVSPBCG 
Its first 3 arguments are the same as for CVSpbcg (see cvspbcg.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSpbcg."""
	ret = cvode.CVBPSpbcg(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpbcg() failed with flag %i"%(ret))
cvode.CVBPSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBPSpbcg.restype = ctypes.c_int

def CVBPSpgmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpgmr links the CVBANDPPRE preconditioner to the CVSPGMR linear solver. It performs the following actions: 
	1) Calls the CVSPGMR specification routine and attaches the CVSPGMR linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPGMR 
	3) Sets the preconditioner setup routine for CVSPGMR 
	4) Sets the preconditioner solve routine for CVSPGMR 
	Its first 3 arguments are the same as for CVSpgmr (see cvspgmr.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSpgmr."""
	ret = cvode.CVBPSpgmr(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpgmr() failed with flag %i"%(ret))
cvode.CVBPSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBPSpgmr.restype = ctypes.c_int

def CVBandPrecFree(bp_data):
 	"""CVBandPrecFree frees the memory allocated by CVBandPrecAlloc in the argument bp_data."""
	ret = cvode.CVBandPrecFree(bp_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecFree() failed with flag %i"%(ret))
cvode.CVBandPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvode.CVBandPrecFree.restype = None

def CVBandPrecGetWorkSpace(bp_data):
	"""CVBandPrecGetWorkSpace returns the CVBandPrec real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVBandPrecGetWorkSpace(bp_data, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVBandPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVBandPrecGetWorkSpace.restype = ctypes.c_int

def CVBandPrecGetNumRhsEvals(bp_data):
 	"""CVBandPrecGetNumGfnEvals returns the number of calls made from CVBANDPRE to the user's RHS."""
	nfevalsBP = ctypes.c_long(0)
	ret = cvode.CVBandPrecGetNumRhsEvals(bp_data, ctypes.byref(nfevalsBP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecGetNumRhsEvals() failed with flag %i"%(ret))
	return nfevalsBP.value
cvode.CVBandPrecGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVBandPrecGetNumRhsEvals.restype = ctypes.c_int

def CVBandPrecGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBANDPRE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVBandPrecGetReturnFlagName(flag)
cvode.CVBandPrecGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVBandPrecGetReturnFlagName.restype = ctypes.c_char_p

################
# cvode_band.h #
################

CVB_MSBJ = 50
CVB_DGMAX = 0.2

class _BandMat(ctypes.Structure):
	_fields_ = [("size", ctypes.c_long), ("mu", ctypes.c_long), ("ml", ctypes.c_long), ("smu", ctypes.c_long), ("data", ctypes.POINTER(ctypes.POINTER(realtype)))]

class _BandMatRow(object):
	"""A class representing a row in a banded matrix"""
	def __init__(self, init, i):
		self.data = init.data
		self.i = i
		self.size = init.size
		self.mu = init.mu
		self.ml = init.ml
		self.smu = init.smu
	
	def __getitem__(self, j):
		"""x.__getitem__(y) <==> x[y]"""
		if (j < 0) or (j >= self.size):
			raise IndexError("Column index out of range")
		y = self.i-j+self.smu
		if ((y < 0) or (y >= self.smu+self.ml+1)):
			return 0.0
		else:
			return self.data.contents.data[j][self.i-j+self.smu]
	
	def __setitem__(self, j, v):
		"""x.__setitem__(i, y) <==> x[i] = y"""
		if (j < 0) or (j >= self.size):
			raise IndexError("Column index out of range")
		y = self.i-j+self.smu
		if ((y < 0) or (y >= self.smu+self.ml+1)):
			if (v != 0):
				raise IndexError("Attempt to store value in unstored element of matrix")
			else:
				return
		self.data.contents.data[j][self.i-j+self.smu] = v

	def __repr__(self):
		ret = "["
		for j in range(self.size):
			y = self.i-j+self.smu
			if ((y < 0) or (y >= self.smu+self.ml+1)):
				ret += "%f"%0.0
			else:
				ret += "%f"%self.data.contents.data[j][y]
			if j != self.size-1:
				ret += ", "
		return ret+"]"

class BandMat(object):
	"""Class representing a banded matrix as a dense matrix"""
	def __init__(self, init, mu = 0, ml = 0, smu = 0):
		"""Instantiates a new dense representation of a banded matrix.
	init	dimension of matrix to create; always square [int]
	mu	width of portion of band above diagonal
	ml	width of portion of band below diagonal
	smu	is the storage upper bandwidth, mu <= smu <= size-1.  The BandGBTRF routine writes the LU factors into the storage for A. The upper triangular factor U, however, may have an upper bandwidth as big as MIN(size-1,mu+ml) because of partial pivoting. The smu field holds the upper bandwidth allocated for A."""
		if type(init) == int:
			self.size = init
			self.mu = mu
			self.ml = ml
			if smu < mu:
				raise ValueError("smu must be greater than or equal to mu")
			self.smu = smu
			self.data = cvode.BandAllocMat(init, mu, ml, smu)
			self.copy = False
		elif type(init) == ctypes.POINTER(_BandMat):
			self.size = init.contents.size
			self.mu = init.contents.mu
			self.ml = init.contents.ml
			self.smu = init.contents.smu
			self.data = init
			self.copy = False
		else:
			raise TypeError("Cannot initialize BandMat from type %s"%(type(init).__name__))
	
	def __getitem__(self, i): #Row major!!!!! index selects a _BandMatRow
		"""x.__getitem__(y) <==> x[y]"""
		if (i < 0) or (i >= self.size):
			raise IndexError("Row index out of range")
		return _BandMatRow(self, i)

	def __setitem__(self, i, row):
		"""x.__setitem__(i, y) <==> x[i] = y"""
		if (i < 0) or (i > self.size):
			raise IndexError("Row index out of range")
		try:
			rowlen = len(row)
		except e:
			raise TypeError("Row assignment requires an object that implements the sequence protocol")
		if rowlen != self.size:
				raise ValueError("Row assignment requires a sequence of length equal to matrix order")
		try:
			r = row[0]
		except:
			raise TypeError("Row assignment requires an object that implements the sequence protocol")
		
		for j in range(self.size):
			y = i-j+self.smu
			if ((y < 0) or (y >= self.smu+self.ml+1)):
				if (row[j] != 0):
					raise IndexError("Attempt to store value in unstored element of matrix")
			self.data.contents.data[j][y] = row[j]
	
	def __repr__(self):
		widest = []
		for j in range(self.size):
			colw = 0
			for i in range(self.mu+self.ml+1):
				colw = max(colw, len("%-10.4f"%(self.data.contents.data[j][i+self.smu])))
			widest.append(colw)

		ret = ""
		for i in range(self.size):
			for j in range(self.size):
				if (i-j+self.smu < 0) or (i-j+self.smu >= self.smu+self.ml+1):
					ret += ("%-"+str(widest[j])+".4f  ")%(0.0)
				else:
					ret += ("%-"+str(widest[j])+".4f  ")%(self.data.contents.data[j][i-j+self.smu])
			ret += "\n"
		return ret

	def __del__(self):
		if not self.copy:
			try:
				cvode.BandFreeMat(self.data)
			except:
				pass

def BandAllocMat(N, mu, ml, smu):
	"""Allocates memory for a banded matrix. Should not be called directly, rather instatiate a BandMat object."""
	return BandMat(N, mu, ml, smu)
cvode.BandAllocMat.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvode.BandAllocMat.restype = ctypes.POINTER(_BandMat)

def BandAllocPiv(N):
	"""BandAllocPiv allocates memory for pivot information to be filled in by the BandGBTRF routine during the factorization of an N by N band matrix. Returns a pointer, which should be passed as is to the [CV]BandPiv* functions.
	N	the size of the banded matrix for which to allocate pivot information space [int]"""
	return cvode.BandAllocPiv(N)
cvode.BandAllocPiv.argtypes = [ctypes.c_int]
cvode.BandAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def BandGBTRF(A, p):
	"""BandGBTRF performs the LU factorization of the N by N band matrix A. This is done using standard Gaussian elimination with partial pivoting.

A successful LU factorization leaves the "matrix" A and the pivot array p with the following information:
	1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.
	2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower triangular matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower triangular part of A contains the multipliers, I-L.

BandGBTRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero.

Important Note: A must be allocated to accommodate the increase in upper bandwidth that occurs during factorization. If mathematically, A is a band matrix with upper bandwidth mu and lower bandwidth ml, then the upper triangular factor U can have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower triangular factor L has lower bandwidth ml. Allocate A with call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are as defined above. The user does not have to zero the "extra" storage allocated for the purpose of factorization. This will handled by the BandGBTRF routine.  ret = cvode.BandGBTRF(A, p)"""
	ret = cvode.BandGBTRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRF() failed with flag %i"%(ret))
cvode.BandGBTRF.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long)]
cvode.BandGBTRF.restype = ctypes.c_long

def BandGBTRS(A, p, b):
	"""BandGBTRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in BandGBTRF. The solution x is returned in b. This routine cannot fail if the corresponding call to BandGBTRF did not fail."""
	ret = cvode.BandGBTRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRS() failed with flag %i"%(ret))
cvode.BandGBTRS.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
cvode.BandGBTRS.restype = None

def BandZero(A):
	"""Zeroes the entire matrix
	A	the matrix [BandMat]"""
	cvode.BandZero(A.data)
cvode.BandZero.argtypes = [_BandMat] 
cvode.BandZero.restype = None

def BandCopy(A, B, copymu, copyml):
	"""Copies the contents of banded matrix B into A, overwriting A's contents.
	A	the destination matrix [BandMat]
	B	the source matrix [BandMat]
	copymu	upper band width of matrices [int]
	copyml	lower band width of matrices [int]"""
	cvode.BandCopy(A.data, B.data, copymu, copyml)
cvode.BandCopy.argtypes = [_BandMat, _BandMat, ctypes.c_long, ctypes.c_long]
cvode.BandCopy.restype = None

def BandScale(c, A):
	"""Scales the matrix A by c
	A	the matrix [BandMat]
	c	the scale factor [float]"""
	cvode.BandScale(c, A.data)
cvode.BandScale.argtypes = [realtype, _BandMat]
cvode.BandScale.restype = None

def BandAddI(A):
	"""Adds 1.0 to the diagonal of the banded matrix A
	A	the matrix [BandMat]"""
	cvode.BandAddI(A.data)
cvode.BandAddI.argtypes = [_BandMat]
cvode.BandAddI.restype = None

def BandFreeMat(A):
	"""BandFreeMat frees the memory allocated by BandAllocMat for the band matrix A.
	A	the matrix [BandMat]"""
	cvode.BandFreeMat(A.data)
	del A.data
cvode.BandFreeMat.argtypes = [_BandMat]
cvode.BandFreeMat.restype = None

def BandFreePiv(p):
	"""BandFreePiv frees the pivot information storage memory p allocated by BandAllocPiv."""
	cvode.BandFreePiv(p)
cvode.BandFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvode.BandFreePiv.restype = None

def BandPrint(A):
	"""Print out the banded matix A"""
	cvode.BandPrint(A.data)
cvode.BandPrint
cvode.BandPrint.restype = None

CVBandJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.POINTER(_BandMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVBandJacFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the Jacobian function. Jacobian functions for banded matrices take N (int = dimension of matrix), muppper (int = upper band width), mlower (int = lower band width), t (float = time step), J (BandMat = Jacobian Matrix), y (NVector), fy (NVector), jac_data (c_void_p), tmp1 (NVector), tmp2 (NVector), and tmp3 (NVector) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVBandJacFn)
	exec 'def __CallbackInterface_%s(N, mupper, mlower, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](N, mupper, mlower, BandMat(J), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVBandJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVBand(cvodememobj, N, mupper, mlower):
	"""A call to the CVBand function links the main CVODE integrator with the CVBAND linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate() 
	N		the size of the ODE system.  
	mupper		the upper bandwidth of the band Jacobian approximation.
	mlower		the lower bandwidth of the band Jacobian approximation."""
	ret = cvode.CVBand(cvodememobj.obj, N, mupper, mlower)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvode.CVBand.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvode.CVBand.restype = ctypes.c_int

def CVBandSetJacFn(cvodememobj, func, jac_data):
	"""CVBandSetJacFn sets the band Jacobian approximation function to be used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func		a python callable suitable for wrapping by WrapCallbackCVBandJacFn
	jac_data	a pointer to user data"""
	ret = cvode.CVBandSetJacFn(cvodememobj.obj, WrapCallbackCVBandJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvode.CVBandSetJacFn.argtypes = [ctypes.c_void_p, CVBandJacFn, ctypes.c_void_p]
cvode.CVBand.restype = ctypes.c_int

def CVBandGetWorkSpace(cvodememobj):
	"""CVBandGetWorkSpace returns the CVBand real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long()
	leniwLS = ctypes.c_long()
	ret = cvode.CVBandGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVBandGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVBandGetWorkSpace.restype = ctypes.c_int

def CVBandGetNumJacEvals(cvodememobj):
 	"""CVBandGetNumGfnEvals returns the number of calls made from CVBAND to the user's jacobian function."""
	retval = ctypes.c_long(0)
	ret = cvode.CVBandGetNumJacEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumJacEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVBandGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVBandGetNumJacEvals.restype = ctypes.c_int

def CVBandGetNumRhsEvals(cvodememobj):
 	"""CVBandGetNumGfnEvals returns the number of calls made from CVBAND to the user's RHS function."""
	retval = ctypes.c_long(0)
	ret = cvode.CVBandGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVBandGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVBandGetNumRhsEvals.restype = ctypes.c_int

def CVBandGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVBAND interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVBandGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVBandGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVBandGetLastFlag.restype = ctypes.c_int

def CVBandGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBAND return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVBandGetReturnFlagName(flag)
cvode.CVBandGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVBandGetReturnFlagName.restype = ctypes.c_char_p

#################
# cvode_dense.h #
#################

class _DenseMat(ctypes.Structure):
	_fields_ = [("M", ctypes.c_long), ("N", ctypes.c_long), ("data", ctypes.POINTER(ctypes.POINTER(realtype)))]

class _DenseMatRow(object):
	"""A class representing a row in a dense matrix"""
	def __init__(self, init, i):
		self.data = init.data
		self.i = i
		self.M = init.M
		self.N = init.N
	
	def __getitem__(self, j):
		"""x.__getitem__(y) <==> x[y]"""
		if (j < 0) or (j >= self.N):
			raise IndexError("Column index out of range")
		else:
			return self.data.contents.data[j][self.i]
	
	def __setitem__(self, j, v):
		"""x.__setitem__(i, y) <==> x[i] = y"""
		if (j < 0) or (j >= self.N):
			raise IndexError("Column index out of range")
		self.data.contents.data[j][self.i] = v

	def __repr__(self):
		ret = "["
		for j in range(self.N):
			ret += "%f"%self.data.contents.data[j][self.i]
			if j != self.N-1:
				ret += ", "
		return ret+"]"

class DenseMat(object):
	"""Class representing a dense matrix"""
	def __init__(self, M, N = None): #M rows by N cols!!!!!
		"""Instantiates a dense matrix object
	M	number of rows
	N	number of columns"""
		if type(M) == int and N is not None and type(N) == int:
			self.M = N
			self.N = M
			self.data = cvode.DenseAllocMat(N,M)
			self.copy = False
		elif type(M) == ctypes.POINTER(_DenseMat):
			self.M = M.contents.M
			self.N = M.contents.N
			self.data = M
			self.copy = True
		else:
			raise TypeError("Cannot initialize DenseMat from type %s"%(type(M).__name__))
	
	def __getitem__(self, index): #row major, i.e. index selects a row
		"""x.__getitem__(y) <==> x[y]; returns row of matrix"""
		if (index >= self.M) or (index < 0):
			raise IndexError
		return _DenseMatRow(self, index)
		
	def __repr__(self):
		widest = []
		for col in range(self.N):
			colw = 0
			for row in range(self.M):
				colw = max(colw, len(str(self.data.contents.data[col][row])))
			widest.append(colw)
		ret = ""
		for row in range(self.M):
			for col in range(self.N):
				ret += ("%-"+str(widest[col])+"f ")%self.data.contents.data[col][row]
			ret += "\n"
		return ret

	def __del__(self):
		if not self.copy:
			try:
				cvode.DenseFreeMat(self.data)
			except:
				pass

def DenseAllocMat(M, N):
	"""Allocates memory for a dense matrix. Should not be called directly, rather instatiate a DenseMat object."""
	return DenseMat(M, N)
cvode.DenseAllocMat.argtypes = [ctypes.c_int, ctypes.c_int]
cvode.DenseAllocMat.restype = ctypes.POINTER(_DenseMat)

def DenseAllocPiv(N):
	"""DenseAllocPiv allocates memory for pivot information to be filled in by the DenseGETRF routine during the factorization of an N by N dense matrix. Returns a pointer, which should be passed as is to the [CV]DensePiv* functions.
	N	the size of the dense matrix for which to allocate pivot information space [int]"""
	return cvode.DenseAllocPiv(N)
cvode.DenseAllocPiv.argtypes = [ctypes.c_int]
cvode.DenseAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def DenseGETRF(A, p):
	"""DenseGETRF performs the LU factorization of the M by N dense matrix A. This is done using standard Gaussian elimination with partial (row) pivoting. Note that this only applies to matrices with M >= N and full column rank.

A successful LU factorization leaves the matrix A and the pivot array p with the following information:
	1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.
	2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower trapezoidal part of A contains the multipliers, I-L.

For square matrices (M=N), L is unit lower triangular.

DenseGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	ret = cvode.DenseGETRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRF() failed with flag %i"%(ret))
cvode.DenseGETRF.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long)]
cvode.DenseGETRF.restype = ctypes.c_long

def DenseGETRS(A, p, b):
	"""DenseGETRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in DenseGETRF. The solution x is returned in b. This routine cannot fail if the corresponding call to DenseGETRF did not fail.
DenseGETRS does NOT check for a squre matrix!"""
	ret = cvode.DenseGETRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRS() failed with flag %i"%(ret))
cvode.DenseGETRS.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
cvode.DenseGETRS.restype = None

def DenseZero(A):
	"""Zeroes the entire matrix
	A	the matrix [DenseMat]"""
	cvode.DenseZero(A.data)
cvode.DenseZero.argtypes = [_DenseMat] 
cvode.DenseZero.restype = None

def DenseCopy(A, B):
	"""Copies the contents of dense matrix B into A, overwriting A's contents.
	A	the destination matrix [DenseMat]
	B	the source matrix [DenseMat]
	copymu	upper band width of matrices [int]
	copyml	lower band width of matrices [int]"""
	cvode.DenseCopy(A.data, B.data)
cvode.DenseCopy.argtypes = [_DenseMat, _DenseMat]
cvode.DenseCopy.restype = None

def DenseScale(c, A):
	"""Scales the matrix A by c
	A	the matrix [DenseMat]
	c	the scale factor [float]"""
	cvode.DenseScale(c, A.data)
cvode.DenseScale.argtypes = [realtype, _DenseMat]
cvode.DenseScale.restype = None

def DenseAddI(A):
	"""Adds 1.0 to the diagonal of the dense matrix A
	A	the matrix [DenseMat]"""
	cvode.DenseAddI(A.data)
cvode.DenseAddI.argtypes = [_DenseMat]
cvode.DenseAddI.restype = None

def DenseFreeMat(A):
	"""DenseFreeMat frees the memory allocated by DenseAllocMat for the band matrix A.
	A	the matrix [DenseMat]"""
	cvode.DenseFreeMat(A.data)
	del A.data
cvode.DenseFreeMat.argtypes = [_DenseMat]
cvode.DenseFreeMat.restype = None

def DenseFreePiv(p):
	"""DenseFreePiv frees the pivot information storage memory p allocated by DenseAllocPiv."""
	cvode.DenseFreePiv(p)
cvode.DenseFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvode.DenseFreePiv.restype = None

def DensePrint(A):
	"""Print out the dense matix A"""
	cvode.DensePrint(A.data)
cvode.DensePrint
cvode.DensePrint.restype = None

CVDenseJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(_DenseMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVDenseJacFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the Jacobian function. Jacobian functions for dense matrices take N (int = dimension of matrix), J (DenseMat = Jacobian Matrix), t (float = time step), y (NVector), fy (NVector), jac_data (c_void_p), tmp1 (NVector), tmp2 (NVector), and tmp3 (NVector) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVDenseJacFn)
	exec 'def __CallbackInterface_%s(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):
	return __ActualCallback[%i](N, DenseMat(J), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVDenseJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVDense(cvodememobj, N):
	"""A call to the CVDense function links the main CVODE integrator with the CVDENSE linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	N		the size of the ODE system."""
	ret = cvode.CVDense(cvodememobj.obj, N)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvode.CVDense.argtypes = [ctypes.c_void_p, ctypes.c_long]
cvode.CVDense.restype = ctypes.c_int

def CVDenseSetJacFn(cvodememobj, func, jac_data):
	"""CVDenseSetJacFn sets the dense Jacobian approximation function to be used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func		a python callable suitable for wrapping by WrapCallbackCVDenseJacFn
	jac_data	a pointer to user data"""
	ret = cvode.CVDenseSetJacFn(cvodememobj.obj, WrapCallbackCVDenseJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvode.CVDenseSetJacFn.argtypes = [ctypes.c_void_p, CVDenseJacFn, ctypes.c_void_p]
cvode.CVDense.restype = ctypes.c_int

def CVDenseGetWorkSpace(cvodememobj):
	"""CVDenseGetWorkSpace returns the CVDense real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long()
	leniwLS = ctypes.c_long()
	ret = cvode.CVDenseGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVDenseGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVDenseGetWorkSpace.restype = ctypes.c_int

def CVDenseGetNumJacEvals(cvodememobj):
 	"""CVDenseGetNumGfnEvals returns the number of calls made from CVDENSE to the user's Jacobian function."""
	retval = ctypes.c_long(0)
	ret = cvode.CVDenseGetNumJacEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumJacEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDenseGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVDenseGetNumJacEvals.restype = ctypes.c_int

def CVDenseGetNumRhsEvals(cvodememobj):
 	"""CVDenseGetNumGfnEvals returns the number of calls made from CVDENSE to the user's RHS function."""
	retval = ctypes.c_long(0)
	ret = cvode.CVDenseGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDenseGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVDenseGetNumRhsEvals.restype = ctypes.c_int

def CVDenseGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVDENSE interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVDenseGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDenseGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVDenseGetLastFlag.restype = ctypes.c_int

def CVDenseGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDENSE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVDenseGetReturnFlagName(flag)
cvode.CVDenseGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVDenseGetReturnFlagName.restype = ctypes.c_char_p

#################
# cvode_spils.h #
#################

CVSPILS_MAXL = 5
CVSPILS_MSBPRE = 50
CVSPILS_DGMAX = 0.2
CVSPILS_DELT = 0.05

CVSPILS_SUCCESS = 0
CVSPILS_MEM_NULL = -1
CVSPILS_LMEM_NULL = -2
CVSPILS_ILL_INPUT = -3
CVSPILS_MEM_FAIL = -4

CVSpilsPrecSetupFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, ctypes.POINTER(ctypes.c_int), realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSetupFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the preconditioner setup function."""
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSetupFn)
	exec 'def __CallbackInterface_%s(t, y, fy, jok, jcurPtr, gamma, P_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(fy), jok, jcurPtr, gamma, P_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSetupFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsPrecSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, realtype, ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSolveFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the preconditioner solve function."""
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSolveFn)
	exec 'def __CallbackInterface_%s(t, y, fy, r, z, gamma, delta, lr, P_data, tmp):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(fy), nvecserial.NVector(r), nvecserial.NVector(z), gamma, delta, lr, P_data, nvecserial.NVector(tmp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsJacTimesVecFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsJacTimesVecFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the jtimes solve function."""
	if func == None:
		return ctypes.cast(None, CVSpilsJacTimesVecFn)
	exec 'def __CallbackInterface_%s(v, Jv, t, y, fy, jac_data, tmp):\n\treturn __ActualCallback[%i](nvecserial.NVector(v), nvecserial.NVector(Jv), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsJacTimesVecFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVSpilsSetPrecType(cvodememobj, pretype):
	"""CVSpilsSetPrecType resets the type of preconditioner, pretype, from the value previously set.  This must be one of PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH."""
	ret = cvode.CVSpilsSetPrecType(cvodememobj.obj, pretype)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPrecType() failed with flag %i"%(ret))
cvode.CVSpilsSetPrecType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVSpilsSetPrecType.restype = ctypes.c_int

def CVSpilsSetGSType(cvodememobj, gstype):
	"""CVSpilsSetGSType specifies the type of Gram-Schmidt orthogonalization to be used. This must be one of the two enumeration constants MODIFIED_GS or CLASSICAL_GS defined in iterative.h. These correspond to using modified Gram-Schmidt and classical Gram-Schmidt, respectively.  Default value is MODIFIED_GS."""
	ret = cvode.CVSpilsSetGSType(cvodememobj.obj, gstype)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetGSType() failed with flag %i"%(ret))
cvode.CVSpilsSetGSType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVSpilsSetGSType.restype = ctypes.c_int

def CVSpilsSetMaxl(cvodememobj, maxl):
	"""CVSpilsSetMaxl resets the maximum Krylov subspace size, maxl, from the value previously set.  An input value <= 0, gives the default value."""
	ret = cvode.CVSpilsSetMaxl(cvodememobj.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetMaxl() failed with flag %i"%(ret))
cvode.CVSpilsSetMaxl.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVSpilsSetMaxl.restype = ctypes.c_int

def CVSpilsSetDelt(cvodememobj, delt):
	"""CVSpilsSetDelt specifies the factor by which the tolerance on the nonlinear iteration is multiplied to get a tolerance on the linear iteration.  Default value is 0.05."""
	ret = cvode.CVSpilsSetDelt(cvodememobj.obj, delt)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetDelt() failed with flag %i"%(ret))
cvode.CVSpilsSetDelt.argtypes = [ctypes.c_void_p, realtype]
cvode.CVSpilsSetDelt.restype = ctypes.c_int

def CVSpilsSetPreconditioner(cvodememobj, pset, psolve, P_data):
	"""CVSpilsSetPreconditioner specifies the PrecSetup and PrecSolve functions.  as well as a pointer to user preconditioner data.  This pointer is passed to PrecSetup and PrecSolve every time these routines are called.  Default is NULL for al three arguments."""
	ret = cvode.CVSpilsSetPreconditioner(cvodememobj.obj, WrapCallbackCVSpilsPrecSetupFn(pset), WrapCallbackCVSpilsPrecSolveFn(psolve), P_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPreconditioner() failed with flag %i"%(ret))
cvode.CVSpilsSetPreconditioner.argtypes = [ctypes.c_void_p, CVSpilsPrecSetupFn, CVSpilsPrecSolveFn, ctypes.c_void_p]
cvode.CVSpilsSetPreconditioner.restype = ctypes.c_int

def CVSpilsSetJacTimesVecFn(cvodememobj, jtimes, jac_data):
	"""CVSpilsSetJacTimesVecFn specifies the jtimes function and a pointer to user Jacobian data. This pointer is passed to jtimes every time the jtimes routine is called.  Default is to use an internal finite difference approximation routine."""
	ret = cvode.CVSpilsSetJacTimesVecFn(cvodememobj.obj, WrapCallbackCVSpilsJacTimesVecFn(jtimes), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetJacTimesVecFn() failed with flag %i"%(ret))
cvode.CVSpilsSetJacTimesVecFn.argtypes = [ctypes.c_void_p, CVSpilsJacTimesVecFn, ctypes.c_void_p]
cvode.CVSpilsSetJacTimesVecFn.restype = ctypes.c_int

def CVSpilsGetWorkSpace(cvodememobj):
	"""CVSpilsGetWorkSpace returns the CVSpils real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVSpilsGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVSpilsGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetWorkSpace.restype = ctypes.c_int

def CVSpilsGetNumPrecEvals(cvodememobj):
	"""CVSpilsGetNumPrecEvals returns the number of preconditioner evaluations, i.e. the number of calls made to PrecSetup with jok==FALSE.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	npevals = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumPrecEvals(cvodememobj.obj, ctypes.byref(npevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecEvals() failed with flag %i"%(ret))
	return npevals.value
cvode.CVSpilsGetNumPrecEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumPrecEvals.restype = ctypes.c_int

def CVSpilsGetNumPrecSolves(cvodememobj):
	"""CVSpilsGetNumPrecSolves returns the number of calls made to PrecSolve.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	npsolves = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumPrecSolves(cvodememobj.obj, ctypes.byref(npsolves))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecSolves() failed with flag %i"%(ret))
	return npsolves.value
cvode.CVSpilsGetNumPrecSolves.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumPrecSolves.restype = ctypes.c_int

def CVSpilsGetNumLinIters(cvodememobj):
	"""CVSpilsGetNumLinIters returns the number of linear iterations.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nliters = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumLinIters(cvodememobj.obj, ctypes.byref(nliters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumLinIters() failed with flag %i"%(ret))
	return nliters.value
cvode.CVSpilsGetNumLinIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumLinIters.restype = ctypes.c_int

def CVSpilsGetNumConvFails(cvodememobj):
	"""CVSpilsGetNumConvFails returns the number of linear convergence failures.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nlcfails = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumConvFails(cvodememobj.obj, ctypes.byref(nlcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumConvFails() failed with flag %i"%(ret))
	return nlcfails.value
cvode.CVSpilsGetNumConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumConvFails.restype = ctypes.c_int

def CVSpilsGetNumJtimesEvals(cvodememobj):
	"""CVSpilsGetNumJtimesEvals returns the number of calls to jtimes.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	njvevals = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumJtimesEvals(cvodememobj.obj, ctypes.byref(njvevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumJtimesEvals() failed with flag %i"%(ret))
	return njvevals.value
cvode.CVSpilsGetNumJtimesEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumJtimesEvals.restype = ctypes.c_int

def CVSpilsGetNumRhsEvals(cvodememobj):
	"""CVSpilsGetNumRhsEvals returns the number of calls to the user f routine due to finite difference Jacobian times vector evaluation.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nfevalsLS = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumRhsEvals(cvodememobj.obj, ctypes.byref(nfevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumRhsEvals() failed with flag %i"%(ret))
	return nfevalsLS.value
cvode.CVSpilsGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumRhsEvals.restype = ctypes.c_int

def CVSpilsGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVSPILS interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	flag = ctypes.c_int(0)
	ret = cvode.CVSpilsGetLastFlag(cvodememobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetLastFlag() failed with flag %i"%(ret))
	return flag.value
cvode.CVSpilsGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVSpilsGetLastFlag.restype = ctypes.c_int

def CVSpilsGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVSPILS return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvode.CVSpilsGetReturnFlagName(flag)
cvode.CVSpilsGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVSpilsGetReturnFlagName.restype = ctypes.c_char_p

##################
# cvode_spbcgs.h #
##################

def CVSpbcg(cvodememobj, pretype, maxl):
	"""A call to the CVSpbcg function links the main CVODE integrator with the CVSPBCG linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in iterative.h. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPBCG solver. Pass 0 to use the default value CVSPBCG_MAXL=5."""
	ret = cvode.CVSpbcg(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpbcg() failed with flag %i"%(ret))
cvode.CVSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSpbcg.restype = ctypes.c_int

#################
# cvode_spgmr.h #
#################

def CVSpgmr(cvodememobj, pretype, maxl):
	"""A call to the CVSpgmr function links the main CVODE integrator with the CVSPGMR linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype	(int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in sundials_iterative.h.  These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPGMR solver. Pass 0 to use the default value CVSPGMR_MAXL=5."""
	ret = cvode.CVSpgmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpgmr() failed with flag %i"%(ret))
cvode.CVSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSpgmr.restype = ctypes.c_int

###################
# cvode_sptfqmr.h #
###################

def CVSptfqmr(cvodememobj, pretype, maxl):
	"""A call to the CVSptfqmr function links the main CVODE integrator with the CVSPTFQMR linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done. This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPTFQMR solver. Pass 0 to use the default value CVSPILS_MAXL=5."""
	ret = cvode.CVSptfqmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSptfqmr() failed with flag %i"%(ret))
cvode.CVSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSptfqmr.restype = ctypes.c_int
