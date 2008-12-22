#cvodes.py is part of the PySUNDIALS package, and is released under the
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

"""Python bindings for the CVODES integrator

The following SUNDIALS header files are wrapped:
	cvodes.h, cvodes_band.h, cvodes_bandpre.h, cvodes_bbdpre.h, cvodes_dense.h,
	cvodes_diag.h, cvodes_spbcgs.h, cvodes_spmgr.h, cvodes_spils.h,
	cvodes_sptfqmr.h
"""

import ctypes
import sundials_core
import nvecserial

realtype = nvecserial.realtype
NVector = nvecserial.NVector

cvodes = sundials_core.loadlib("cvodes")

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
CV_EE = 4

#Internal task job type
CV_NORMAL = 1
CV_ONE_STEP = 2
CV_NORMAL_TSTOP = 3
CV_ONE_STE_TSTOP = 4

#ism
CV_SIMULTANEOUS = 1
CV_STAGGERED = 2
CV_STAGGERED1 = 3

#DQtype
CV_CENTERED = 1
CV_FORWARD = 2

#interp
CV_HERMITE = 1
CV_POLYNOMIAL = 2

#CVodes() return flags
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

CV_NO_QUAD = -30
CV_QRHSFUNC_FAIL = -31
CV_FIRST_QRHSFUNC_ERR = -32
CV_REPTD_QRHSFUNC_ERR = -33
CV_UNREC_QRHSFUNC_ERR = -34

CV_BAD_IS = -40
CV_NO_SENS = -41
CV_SRHSFUNC_FAIL = -42
CV_FIRST_SRHSFUNC_ERR = -43
CV_REPTD_SRHSFUNC_ERR = -44
CV_UNREC_SRHSFUNC_ERR = -45

#CVODEA return flags

CV_ADJMEM_NULL = -101
CV_BAD_TB0 = -103
CV_BCKMEM_NULL = -104
CV_REIFWD_FAIL = -105
CV_FWD_FAIL = -106
CV_BAD_ITASK = -107
CV_BAD_TBOUT = -108
CV_GETY_BADT = -109

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
			cvodes.CVodeFree(ctypes.byref(p))
cvodes.CVodeFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvodes.CVodeFree.restype = None

CVRhsFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVRhsFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a RHS function.

The callable python object must take exactly 4 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
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
	t (realtype)		the current value of the independent variable
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

The callable python object must take exactly 5 parameters, which will be passed in as
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

CVQuadRhsFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVQuadRhsFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a Quaderature RHS function.

The callable python object must take exactly 4 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
	y (NVector)		the vector of current dependent values
	yQdot (NVector)		undefined values, contents should be set to the new values of y
	fQ_data (c_void_p)	pointer to user data set by CVodeSetQuadFdata

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVQuadRhsFn)
	exec 'def __CallbackInterface_%s(t, y, yQdot, fQ_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yQdot), fQ_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVQuadRhsFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSensRhsFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSensRhsFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly CVodesetSensRhsFn

The callable python object must take exactly 6 parameters, which will be passed in as
	Ns (int)		the number of sensitivities
	t (realtype)		the current value of the independent variable
	y (NVector)		the vector of current dependent values
	yS (NVector)		the vector of current sensitivity values
	ySdot (NVector)		undefined values, contents should be set to the new values of yS
	fS_data (c_void_p)	pointer to user data set by CVodeSetSensFdata

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVSensRhsFn)
	exec 'def __CallbackInterface_%s(Ns, t, y, ydot, yS, ySdot, fS_data, tmp1, tmp2):\n\treturn __ActualCallback[%i](Ns, t, nvecserial.NVector(y), nvecserial.NVector(ydot), nvecserial.NVector(yS), nvecserial.NVector(ySdot), fS_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSensRhsFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSensRhs1Fn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSensRhs1Fn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly CVodesetSensRhsFn

The callable python object must take exactly 10 parameters, which will be passed in as
	Ns (int)		the number of sensitivities
	t (realtype)		the current value of the independent variable
	y (NVector)		the vector of current dependent values
	ydot (NVector)		undefined values, contents should be set to the new values of y
	iS (int)		the current sensitivity
	yS (NVector)		the vector of current sensitivity values
	ySdot (NVector)		undefined values, contents should be set to the new values of yS
	fS_data (c_void_p)	pointer to user data set by CVodeSetSensFdata
	tmp1 (NVector)		preallocated temporary storage
	tmp2 (NVector)		preallocated temporary storage

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVSensRhs1Fn)
	exec 'def __CallbackInterface_%s(Ns, t, y, ydot, iS, yS, ySdot, fS_data, tmp1, tmp2):\n\treturn __ActualCallback[%i](Ns, t, nvecserial.NVector(y), nvecserial.NVector(ydot), iS, nvecserial.NVector(yS), nvecserial.NVector(ySdot), fS_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSensRhs1Fn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVRhsFnB = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVRhsFnB(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a RHS function.

The callable python object must take exactly 5 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
	y (NVector)		the vector of current dependent values
	yB (NVector)		the vector of current dependent values
	yBdot (NVector)		undefined values, contents should be set to the new values of y
	f_data (c_void_p)	pointer to user data set by CVodeSetFdata

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVRhsFnB)
	exec 'def __CallbackInterface_%s(t, y, yB, yBdot, f_dataB):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(yBdot), f_dataB)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVRhsFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVQuadRhsFnB = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVQuadRhsFnB(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by any functions that specify a Quaderature RHS function.

The callable python object must take exactly 5 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
	y (NVector)		the vector of current dependent values
	yB (NVector)		the vector of current dependent values
	qBdot (NVector)		undefined values, contents should be set to the new values of y
	fQ_dataB (c_void_p)	pointer to user data set by CVodeSetQuadFdata

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVQuadRhsFnB)
	exec 'def __CallbackInterface_%s(t, y, yB, qBdot, fQ_dataB):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(qBdot), fQ_dataB)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVQuadRhsFnB(eval("__CallbackInterface_%s"%(func.func_name)))
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
	obj = cvodes.CVodeCreate(lmm, iter)
	if obj == None:
		raise AssertionError("SUNDIALS ERROR: CVodeCreate() failed - returned NULL pointer")
	return CVodeMemObj(obj)
cvodes.CVodeCreate.argtypes = [ctypes.c_int, ctypes.c_int]
cvodes.CVodeCreate.restype = ctypes.c_void_p

def CVodeSetErrHandlerFn(cvodememobj, func,  eh_data):
	"""Sets a user provided error handling function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func	is a python callable taking the parameters
		error_code (int)	the code of the error
		module (string)		the name of the module producing the error
		function_name (string)	the function name in which the error appeared
		message (string)	a message describing the nature of the error
		eh_data (c_void_p)	a pointer to user data as set by eh_data
	
	eh_data (c_void_p)	a pointer to any user data you would like passed through to your custom handler. The specified function will be passed eh_data as its eh_data parameter."""
	ret = cvodes.CVodeSetErrHandlerFn(cvodememobj.obj, WrapCallbackCVErrHandlerFn(func), eh_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrHanderFn() failed with flag %i"%(ret))
cvodes.CVodeSetErrHandlerFn.argtypes = [ctypes.c_void_p, CVErrHandlerFn, ctypes.c_void_p]
cvodes.CVodeSetErrHandlerFn.restype = ctypes.c_int

def CVodeSetErrFile(cvodememobj, fileobj):
	"""Sets the file where all error messages and warnings will be written if the default error handling function is used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	fileobj (file object)		must be open for writing"""
	ret = cvodes.CVodeSetErrFile(cvodememobj.obj, sundials_core.fdopen(fileobj.fileno, fileobj.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrFile() failed with flag %i"%(ret))
cvodes.CVodeSetErrFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetErrFile.restype = ctypes.c_int

def CVodeSetFdata(cvodememobj, f_data):
	"""Sets the pointer to any user data. If called, whenever the user's right hand side function is evaluated, it will be passed f_data as its f_data parameter."""
	ret = cvodes.CVodeSetFdata(cvodememobj.obj, f_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetFdata() failed with flag %i"%(ret))
cvodes.CVodeSetFdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetFdata.restype = ctypes.c_int

def CVodeSetEwtFn(cvodememobj, func, e_data):
	"""Sets a user provided error weight function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func	is a python callable taking the following parameters which must return 0 in the case of success, otherwise non-zero
		y (NVector)		the vector of current dependent variables
		ewt (NVector)		undefined values, contents should be set to the error weights
		e_data (c_void_p)	a pointer to user data as set by e_data

	e_data (c_void_p)	is a pointer to any user data you would like passed through to your custom function. The specified function will be passed e_data as its e_data parameter."""
	ret = cvodes.CVodeSetEwtFn(cvodememobj.obj, WrapCallbackCVEwtFn(func), e_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetEwtFn() failed with flag %i"%(ret))
cvodes.CVodeSetEwtFn.argtypes = [ctypes.c_void_p, CVEwtFn, ctypes.c_void_p]
cvodes.CVodeSetEwtFn.restype = ctypes.c_int

def CVodeSetMaxOrd(cvodememobj, maxord):
	"""Sets the maximum lmm order to be used by the solver. Defaults to 12 for CV_ADAMS, and 5 for CV_BDF
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxord (int)"""
	ret = cvodes.CVodeSetMaxOrd(cvodememobj.obj, maxord)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxOrd() failed with flag %i"%(ret))
cvodes.CVodeSetMaxOrd.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxOrd.restype = ctypes.c_int

def CVodeSetMaxNumSteps(cvodememobj, maxsteps):
	"""Sets the maximum number of internal steps to be taken by the solver in its attempt to reach tout.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxsteps (int)"""
	ret = cvodes.CVodeSetMaxNumSteps(cvodememobj.obj, maxsteps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNumSteps() failed with flag %i"%(ret))
cvodes.CVodeSetMaxNumSteps.argtypes = [ctypes.c_void_p, ctypes.c_long] 
cvodes.CVodeSetMaxNumSteps.restype = ctypes.c_int

def CVodeSetMaxHnilWarns(cvodememobj, mxhnil):
	"""Sets the maximum number of warning messages issued by the solver that t+h==t on the next internal step. A value of -1 means no such warnings are issued. Default is 10.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	mxhnil (int)"""
	ret = cvodes.CVodeSetMaxHnilWarns(cvodememobj.obj, mxhnil)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxHnilWarns() failed with flag %i"%(ret))
cvodes.CVodeSetMaxHnilWarns.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxHnilWarns.restype = ctypes.c_int

def CVodeSetStabLimDet(cvodememobj, stldet):
	"""Turns stability limit detection on or off
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	stldet (int)	0 = OFF, 1 = ON

When CV_BDF is used and order is 3 or greater, CVsldet is called to detect stability limit. If limit is detected, the order is reduced. Default is off"""
	ret = cvodes.CVodeSetStabLimDet(cvodememobj.obj, stldet)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStabLimDet() failed with flag %i"%(ret))
cvodes.CVodeSetStabLimDet.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetStabLimDet.restype = ctypes.c_int

def CVodeSetInitStep(cvodememobj, hin):
	"""Sets initial step size. CVODE estimates this value by default.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	hin (float)"""
	ret = cvodes.CVodeSetInitStep(cvodememobj.obj, hin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetInitStep() failed with flag %i"%(ret))
cvodes.CVodeSetInitStep.argtypes = [ctypes.c_void_p, realtype] 
cvodes.CVodeSetInitStep.restype = ctypes.c_int

def CVodeSetMinStep(cvodememobj, hmin):
	"""Sets the absolute minimum of step size allowed. Default is 0.0
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	hmin (float)"""
	ret = cvodes.CVodeSetMinStep(cvodememobj.obj, hmin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMinStep() failed with flag %i"%(ret))
cvodes.CVodeSetMinStep.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetMinStep.restype = ctypes.c_int

def CVodeSetMaxStep(cvodememobj, hmax):
	"""Sets the maximum absolute value of step size allowed. Default is infinity.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	hmax (float)"""
	ret = cvodes.CVodeSetMaxStep(cvodememobj.obj, hmax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxStep() failed with flag %i"%(ret))
cvodes.CVodeSetMaxStep.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetMaxStep.restype = ctypes.c_int

def CVodeSetStopTime(cvodememobj, tstop):
	"""Sets the independent variable value past which the solution is not to proceed. Default is infinity.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	tstop (float)"""
	ret = cvodes.CVodeSetStopTime(cvodememobj.obj, tstop)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStopTime() failed with flag %i"%(ret))
cvodes.CVodeSetStopTime.argtypes = [ctypes.c_void_p, realtype] 
cvodes.CVodeSetStopTime.restype = ctypes.c_int

def CVodeSetMaxErrTestFails(cvodememobj, maxnef):
	"""Sets the maximum number of error test failures in attempting one step. Default is 7.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxnef (int)"""
	ret = cvodes.CVodeSetMaxErrTestFails(cvodememobj.obj, maxnef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxErrTestFails() failed with flag %i"%(ret))
cvodes.CVodeSetMaxErrTestFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxErrTestFails.restype = ctypes.c_int

def CVodeSetMaxNonlinIters(cvodememobj, maxcor):
	"""Sets the maximum number of nonlinear solver iterations at one solution. Default is 3.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxcor (int)"""
	ret = cvodes.CVodeSetMaxNonlinIters(cvodememobj.obj, maxcor)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNonlinIters() failed with flag %i"%(ret))
cvodes.CVodeSetMaxNonlinIters.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxNonlinIters.restype = ctypes.c_int

def CVodeSetMaxConvFails(cvodememobj, maxncf):
	"""Sets the maximum number of convergence failures allowed in attempting one step. Default is 10.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxncf (int)"""
	ret = cvodes.CVodeSetMaxConvFails(cvodememobj.obj, maxncf)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxConvFails() failed with flag %i"%(ret))
cvodes.CVodeSetMaxConvFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxConvFails.restype = ctypes.c_int

def CVodeSetNonlinConvCoef(cvodememobj, nlscoef):
	"""Sets the coefficient in the nonlinear convergence test. Default is 0.1
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	nlscoef (float)"""
	ret = cvodes.CVodeSetNonlinConvCoef(cvodememobj.obj, nlscoef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetNonlinConvCoef() failed with flag %i"%(ret))
cvodes.CVodeSetNonlinConvCoef.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetNonlinConvCoef.restype = ctypes.c_int

def CVodeSetIterType(cvodememobj, iter):
	"""Changes the current nonlinear iteration type. Initially set by CVodeCreate().
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	iter (int)	specifies whether functional or newton iteration will be used. Legal values are: 
		CV_FUNCTIONAL	does not require linear algebra 
		CV_NEWTON	requires the solution of linear systems. Requires the specification of a CVODE linear solver. (Recommended for stiff problems)."""
	ret = cvodes.CVodeSetIterType(cvodememobj.obj, iter)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetIterType() failed with flag %i"%(ret))
cvodes.CVodeSetIterType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetIterType.restype = ctypes.c_int

def CVodeSetTolerances(cvodememobj, itol, reltol, abstol):
	"""Changes the integration tolerances between calls to CVode(). Initially (re)set by CVodeMalloc()/CVodeReInit().
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
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
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvodes.CVodeSetTolerances(cvodememobj.obj, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS, CV_SV or CV_WF")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetTolerances() failed with flag %i"%(ret))
cvodes.CVodeSetTolerances.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeSetTolerances.restype = ctypes.c_int

def CVodeMalloc(cvodememobj, func, t0, y0, itol, reltol, abstol):
	"""CVodeMalloc allocates and initializes memory for a problem to be solved by CVODE.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	func (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		t (realtype)		the current value of the independent variable
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
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvodes.CVodeMalloc(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeMalloc() failed with flag %i"%(ret))
	cvodememobj.dealloc = True
cvodes.CVodeMalloc.argtypes = [ctypes.c_void_p, CVRhsFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeMalloc.restype = ctypes.c_int

def CVodeReInit(cvodememobj, func, t0, y0, itol, reltol, abstol):
	"""CVodeReInit re-initializes CVode for the solution of a problem, where a prior call to CVodeMalloc has been made with the same problem size N. CVodeReInit performs the same input checking and initializations that CVodeMalloc does.  But it does no memory allocation, assuming that the existing internal memory is sufficient for the new problem.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	func (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		t (realtype)		the current value of the independent variable
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
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(abstol))
		elif type(abstol) == float or type(abstol) == int:
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, ctypes.byref(realtype(abstol)))
		elif abstol == None:
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be a floating point number if itol is CV_SS")
	elif itol == CV_SV:
		if type(abstol) == NVector:
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol.data)
		elif abstol == None:
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be an NVector if itol is CV_SV")
	elif itol == CV_WF:
		if abstol == None:
			ret = cvodes.CVodeReInit(cvodememobj.obj, WrapCallbackCVRhsFn(func), t0, y0.data, itol, reltol, abstol)
		else:
			raise TypeError("abstol must be None if itol is CV_WF")
	else:
		raise ValueError("itol must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeReInit() failed with flag %i"%(ret))
cvodes.CVodeReInit.argtypes = [ctypes.c_void_p, CVRhsFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeReInit.restype = ctypes.c_int

def CVodeRootInit(cvodememobj, nrtfn, func, g_data):
	"""CVodeRootInit initializes a rootfinding problem to be solved during the integration of the ODE system. It must be called after CVodeCreate, and before CVode. The arguments are:
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	nrtfn (int)			number of functions whose roots are to be found. nrtfn must be positive.
	func (function)			a python callable which will be automatically wrapped by WrapCallbackCVRootFn() which defines the functions g_i (0 < i < nrtfn) whose roots are sought. The function should take the the following parameters
		t (realtype)		the current value of the independent variable
		y (NVector)		the vector of current dependent values
		gout (NVector)		undefined values, contents should be set to the roots desired. If any element of gout is zero upon return, CVode will return early specifying the numbers of roots found, i.e. the number of 0 elements of gout
		g_data (c_void_p)	pointer to user data set by CVodeRootInit
	g_data (c_void_p)		a pointer to user data that will be passed to the user's root evaluation function every time it is called. The function receives the value of g_data in its g_data parameter

If a new problem is to be solved with a call to CVodeReInit, where the new problem has no root functions but the prior one did, then call CVodeRootInit with nrtfn = 0."""
	ret = cvodes.CVodeRootInit(cvodememobj.obj, nrtfn, WrapCallbackCVRootFn(func), g_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeRootInit() failed with flag %i"%(ret))
cvodes.CVodeRootInit.argtypes = [ctypes.c_void_p, ctypes.c_int, CVRootFn, ctypes.c_void_p]
cvodes.CVodeRootInit.restype = ctypes.c_int

def CVodeSetQuadFdata(cvodememobj, fQ_data):
	"""Sets the pointer to any user data. If called, whenever the user's right hand side function is evaluated, it will be passed fQ_data as its fQ_data parameter."""
	ret = cvodes.CVodeSetQuadFdata(cvodememobj.obj, fQ_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetQuadFdata() failed with flag %i"%(ret))
cvodes.CVodeSetQuadFdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetQuadFdata.restype = ctypes.c_int

def CVodeSetQuadErrCon(cvodememobj, errconQ, itolQ, reltolQ, abstolQ):
	"""Specify the consideration of quadrature variables in error control
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	errconQ (int)			1 or 0, 1 indicates quadrature cvariables are considered
	itolQ (int)			one of CV_SS or CV_SV
	reltolQ (float)			the relative tolerance scalar
	abstolQ (float/NVector)		absolute tolerance(s)"""
	if itolQ == CV_SS:
		if type(abstolQ) == realtype:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, ctypes.byref(abstolQ))
		elif type(abstolQ) == float or type(abstolQ) == int:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, ctypes.byref(realtype(abstolQ)))
		elif abstolQ == None:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, abstolQ)
		else:
			raise TypeError("abstolQ must be a floating point number if itolQ is CV_SS")
	elif itolQ == CV_SV:
		if type(abstolQ) == NVector:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, abstolQ.data)
		elif abstolQ == None:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, abstolQ)
		else:
			raise TypeError("abstolQ must be an NVector if itol is CV_SV")
	elif itolQ == CV_WF:
		if abstolQ == None:
			ret = cvodes.CVodeSetQuadErrCon(cvodememobj.obj, errconQ, itolQ, reltolQ, abstolQ)
		else:
			raise TypeError("abstolQ must be None if itolQ is CV_WF")
	else:
		raise ValueError("itolQ must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetQuadErrCon() failed with flag %i"%(ret))
cvodes.CVodeSetQuadErrCon.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeSetQuadErrCon.restype = ctypes.c_int

def CVodeQuadMalloc(cvodememobj, fQ, yQ0):
	"""CVodeMalloc allocates and initializes memory for a problem to be solved by CVODE.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	fQ (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		t (realtype)		the current value of the independent variable
		y (NVector)		the vector of current dependent values
		yQdot (NVector)		undefined values, contents should be set to the new values of y
		fQ_data (c_void_p)	pointer to user data set by CVodeSetQuadFdata
	yQ0 (NVector)			the initial quadrature vector yQ(t0)."""
	ret = cvodes.CVodeQuadMalloc(cvodememobj.obj, WrapCallbackCVQuadRhsFn(fQ), yQ0.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeQuadMalloc() failed with flag %i"%(ret))
cvodes.CVodeQuadMalloc.argtypes = [ctypes.c_void_p, CVQuadRhsFn, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeQuadMalloc.restype = ctypes.c_int

def CVodeQuadReInit(cvodememobj, fQ, yQ0):
	"""CVodeQuadReInit reinitialises quadrature values.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	fQ (function)			a python callable defining the right hand side function in y' = f(t,y), which is automatically wrapped by WrapCallbackCVRhsFn. The function specified takes the following parameters
		t (realtype)		the current value of the independent variable
		y (NVector)		the vector of current dependent values
		yQdot (NVector)		undefined values, contents should be set to the new values of y
		fQ_data (c_void_p)	pointer to user data set by CVodeSetQuadFdata
	yQ0 (NVector)			the initial quadrature vector yQ(t0)."""
	ret = cvodes.CVodeQuadReInit(cvodememobj.obj, WrapCallbackCVQuadRhsFn(fQ), yQ0.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeQuadReInit() failed with flag %i"%(ret))
cvodes.CVodeQuadReInit.argtypes = [ctypes.c_void_p, CVQuadRhsFn, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeQuadReInit.restype = ctypes.c_int

def CVodeSetSensRhsFn(cvodememobj, f, fS_dataS):
	ret = cvodes.CVodeSetSensRhsFn(cvodememobj.obj, WrapCallbackCVSensRhsFn(f), fS_dataS)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensRhsFn() failed with flag %i"%(ret))
cvodes.CVodeSetSensRhsFn.argtypes = [ctypes.c_void_p, CVSensRhsFn, ctypes.c_void_p]
cvodes.CVodeSetSensRhsFn.restype = ctypes.c_int

def CVodeSetSensRhs1Fn(cvodememobj, fS, fS_data):
	ret = cvodes.CVodeSetSensRhs1Fn(cvodememobj.obj, WrapCallbackCVSensRhs1Fn(fS), fS_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensRhs1Fn() failed with flag %i"%(ret))
cvodes.CVodeSetSensRhs1Fn.argtypes = [ctypes.c_void_p, CVSensRhs1Fn, ctypes.c_void_p]
cvodes.CVodeSetSensRhs1Fn.restype = ctypes.c_int

def CVodeSetSensDQMethod(cvodememobj, DQtype, DQrhomax):
	"""Controls the selection of finite difference schemses used in evaluating the sensitivity right hand sides.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	DQtype (int)		one of CV_CENETERED (default), CV_FORWARD, CV_SIMULTANEOUS or CV_SEPARATE
	DQrhormax (float)	default 0"""
	ret = cvodes.CVodeSetSensDQMethod(cvodememobj.obj, DQtype, DQrhomax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensDQMethod() failed with flag %i"%(ret))
cvodes.CVodeSetSensDQMethod.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype]
cvodes.CVodeSetSensDQMethod.restype = ctypes.c_int

def CVodeSetSensErrCon(cvodememobj, errconS):
	"""Consider sensitivity variables in error control?
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	errconS (int)				1 = true"""
	ret = cvodes.CVodeSetSensErrCon(cvodememobj.obj, errconS)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensErrCon() failed with flag %i"%(ret))
cvodes.CVodeSetSensErrCon.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetSensErrCon.restype = ctypes.c_int

def CVodeSetSensMaxNonlinIters(cvodememobj, maxcorS):
	"""Sets the maximum number of nonlinear solver iterations at one solution. Default is 3.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	maxcor (int)"""
	ret = cvodes.CVodeSetSensMaxNonlinIters(cvodememobj.obj, maxcorS)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensMaxNonlinIters() failed with flag %i"%(ret))
cvodes.CVodeSetSensMaxNonlinIters.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetSensMaxNonlinIters.restype = ctypes.c_int

def CVodeSetSensParams(cvodememobj, p, pbar, plist):
	"""Set parameter information
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	p (*realtype)				pointer to problem paramters (usually part of some user data structure)
	pbar (*realtype)			
	plist (*realtype)			
	if pbar is not None:
		wpbar = (realtype*len(pbar))()
		for i in range(len(pbar)):
			wpbar[i] = pbar[i]
	else:
		wpbar = pbar
	if plist is not None:
		wplist = (ctypes.c_int*len(plist))()
		for i in range(len(plist)):
			wplist[i] = plist[i]
	else:
		wplist = plist
	ret = cvodes.CVodeSetSensParams(cvodememobj.obj, p, wpbar, wplist)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensParams() failed with flag %i"%(ret))
cvodes.CVodeSetSensParams.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype), ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int)]
cvodes.CVodeSetSensParams.restype = ctypes.c_int

def CVodeSetSensTolerances(cvodememobj, itolS, reltolS, abstolS):
	if itolS == CV_SS:
		if type(abstolS) == realtype:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, ctypes.byref(abstolS))
		elif type(abstolS) == float or type(abstolS) == int:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, ctypes.byref(realtype(abstolS)))
		elif abstolS == None:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, abstolS)
		else:
			raise TypeError("abstolS must be a floating point number if itolS is CV_SS")
	elif itolS == CV_SV:
		if type(abstolS) == NVector:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, abstolS.data)
		elif abstolS == None:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, abstolS)
		else:
			raise TypeError("abstolS must be an NVector if itolS is CV_SV")
	elif itolS == CV_WF:
		if abstolS == None:
			ret = cvodes.CVodeSetSensTolerances(cvodememobj.obj, itolS, reltolS, abstolS)
		else:
			raise TypeError("abstolS must be None if itolS is CV_WF")
	else:
		raise ValueError("itolS must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetSensTolerances() failed with flag %i"%(ret))
cvodes.CVodeSetSensTolerances.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeSetSensTolerances.restype = ctypes.c_int

def CVodeSensMalloc(cvodememobj, Ns, ism, yS0):
	ret = cvodes.CVodeSensMalloc(cvodememobj.obj, Ns, ism, yS0.cdata)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSensMalloc() failed with flag %i"%(ret))
cvodes.CVodeSensMalloc.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(nvecserial._NVector))]
cvodes.CVodeSensMalloc.restype = ctypes.c_int

def CVodeSensReInit(cvodememobj, ism, yS0):
	ret = cvodes.CVodeSensReInit(cvodememobj.obj, ism, yS0.cdata)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSensReInit() failed with flag %i"%(ret))
cvodes.CVodeSensReInit.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(nvecserial._NVector))]
cvodes.CVodeSensReInit.restype = ctypes.c_int

def CVodeSensToggleOff(cvodememobj):
	ret = cvodes.CVodeSensToggleOff(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSensToggleOff() failed with flag %i"%(ret))
cvodes.CVodeSensToggleOff.argtypes = [ctypes.c_void_p]
cvodes.CVodeSensToggleOff.restype = ctypes.c_int

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
	ret = cvodes.CVode(cvodememobj.obj, tout, yout.data, tret, itask)
	return ret
cvodes.CVode.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype), ctypes.c_int]
cvodes.CVode.restype = ctypes.c_int

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
	ret = cvodes.CVodeGetDky(cvodememobj.obj, t, k, dky.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetDky() failed with flag %i"%(ret))
cvodes.CVodeGetDky.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)] 
cvodes.CVodeGetDky.restype = ctypes.c_int

def CVodeGetWorkSpace(cvodememobj):
	"""CVodeGetWorkSpace returns the CVODE real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrw = ctypes.c_long()
	leniw = ctypes.c_long()
	ret = cvodes.CVodeGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrw), ctypes.byref(leniw))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetWorkSpace() failed with flag %i"%(ret))
	return (lenrw, leniw)
cvodes.CVodeGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetWorkSpace.restype = ctypes.c_int

def CVodeGetNumSteps(cvodememobj):
	"""CVodeGetNumSteps returns the cumulative number of internal steps taken by the solver
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvodes.CVodeGetNumSteps(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSteps() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumSteps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSteps.restype = ctypes.c_int

def CVodeGetNumRhsEvals(cvodememobj):
	"""CVodeGetNumRhsEvals returns the number of calls to the user's f function
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvodes.CVodeGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumRhsEvals.restype = ctypes.c_int

def CVodeGetNumLinSolvSetups(cvodememobj):
	"""CVodeGetNumLinSolvSetups returns the number of calls made to the linear solver's setup routine
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvodes.CVodeGetNumLinSolvSetups(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumLinSolvSetups() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumLinSolvSetups.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumLinSolvSetups.restype = ctypes.c_int

def CVodeGetNumErrTestFails(cvodememobj):
	"""CVodeGetNumErrTestFails returns the number of local error test failures that have occured
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvodes.CVodeGetNumErrTestFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumErrTestFails() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumErrTestFails.restype = ctypes.c_int

def CVodeGetLastOrder(cvodememobj):
	"""CVodeGetLastOrder returns the order used during the last internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvodes.CVodeGetLastOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastOrder() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetLastOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVodeGetLastOrder.restype = ctypes.c_int

def CVodeGetCurrentOrder(cvodememobj):
	"""CVodeGetCurrentOrder returns the order to be used on the next internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvodes.CVodeGetCurrentOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentOrder() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetCurrentOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVodeGetCurrentOrder.restype = ctypes.c_int

def CVodeGetNumStabLimOrderReds(cvodememobj):
	"""CVodeGetNumStabLimOrderReds returns the number of order reductions due to stability limit detection
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvodes.CVodeGetNumStabLimOrderReds(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumStabLimOrderReds() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumStabLimOrderReds.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumStabLimOrderReds.restype = ctypes.c_int

def CVodeGetActualInitStep(cvodememobj):
	"""CVodeGetActualInitStep returns the actual initial step size used by CVODE
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvodes.CVodeGetActualInitStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetActualInitStep() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetActualInitStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvodes.CVodeGetActualInitStep.restype = ctypes.c_int

def CVodeGetLastStep(cvodememobj):
	"""CVodeGetLastStep returns the step size for the last internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvodes.CVodeGetLastStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastStep() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetLastStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvodes.CVodeGetLastStep.restype = ctypes.c_int

def CVodeGetCurrentStep(cvodememobj):
	"""CVodeGetCurrentStep returns the step size to be attempted on the next internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvodes.CVodeGetCurrentStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentStep() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetCurrentStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvodes.CVodeGetCurrentStep.restype = ctypes.c_int

def CVodeGetCurrentTime(cvodememobj):
	"""CVodeGetCurrentTime returns the current internal time reached by the solver
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvodes.CVodeGetCurrentTime(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentTime() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetCurrentTime.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvodes.CVodeGetCurrentTime.restype = ctypes.c_int

def CVodeGetTolScaleFactor(cvodememobj):
	"""CVodeGetTolScaleFactor returns a suggested factor by which the user's tolerances should be scaled when too much accuracy has been requested for some internal step
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvodes.CVodeGetTolScaleFactor(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetTolScaleFactor() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetTolScaleFactor.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvodes.CVodeGetTolScaleFactor.restype = ctypes.c_int

def CVodeGetErrWeights(cvodememobj, eweight):
	"""CVodeGetErrWeights returns the current error weight vector in eweight.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	eweight	(NVector)		is set to the current error weight vector"""
	if type(cvodememobj).__name__ == 'CVodeMemObj':
		cvo = cvodememobj.obj
	else:
		cvo = cvodememobj
	ret = cvodes.CVodeGetErrWeights(cvo, eweight.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetErrWeights() failed with flag %i"%(ret))
cvodes.CVodeGetErrWeights.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetErrWeights.restype = ctypes.c_int

def CVodeGetEstLocalErrors(cvodememobj, ele):
	"""CVodeGetEstLocalErrors returns the vector of estimated local errors in ele.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	ele (NVector)			is set to the vector of current estimated local errors"""
	ret = cvodes.CVodeGetEstLocalErrors(cvodememobj.obj, ele.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetEstLocalErrors() failed with flag %i"%(ret))
cvodes.CVodeGetEstLocalErrors.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetEstLocalErrors.restype = ctypes.c_int

def CVodeGetNumGEvals(cvodememobj):
	"""CVodeGetNumGEvals returns the number of calls to the user's g function (for rootfinding)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumGEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumGEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumGEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumGEvals.restype = ctypes.c_int

def CVodeGetRootInfo(cvodememobj, numroots):
	"""CVodeGetRootInfo returns the indices for which the function g_i was found to have a root. A list of zeroes and ones is returned, where a one in index postion i indicates a root for g_i was found.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	numroots (int)			the number of functions for which roots should be found."""
	rootsfound = (ctypes.c_int*numroots)()
	ret = cvodes.CVodeGetRootInfo(cvodememobj.obj, rootsfound)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVode() failed with flag %i"%(ret))
	return [rootsfound[i] for i in range(numroots)]
cvodes.CVodeGetRootInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVodeGetRootInfo.restype = ctypes.c_int
	
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
	ret = cvodes.CVodeGetIntegratorStats(cvodememobj.obj, ctypes.byref(nsteps), ctypes.byref(nfevals), ctypes.byref(nlinsetups), ctypes.byref(netfails), ctypes.byref(qlast), ctypes.byref(qcur), ctypes.byref(hinused), ctypes.byref(hlast), ctypes.byref(hcur), ctypes.byref(tcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetIntegratorStats() failed with flag %i"%(ret))
	return (nsteps.value, nfevals.value, nlinsetups.value, netfails.value, qlast.value, qcur.value, hinused.value, hlast.value, hcur.value, tcur.value)
cvodes.CVodeGetIntegratorStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(realtype), ctypes.POINTER(realtype), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
cvodes.CVodeGetIntegratorStats.restype = ctypes.c_int

def CVodeGetNumNonlinSolvIters(cvodememobj):
	"""CVodeGetNumNonlinSolvIters returns the number of nonlinear solver iterations performed.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumNonlinSolvIters(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvIters() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumNonlinSolvIters.restype = ctypes.c_int

def CVodeGetNumNonlinSolvConvFails(cvodememobj):
	"""CVodeGetNumNonlinSolvConvFails returns the number of nonlinear convergence failures.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumNonlinSolvConvFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvConvFails() failed with flag %i"%(ret))
	return retval.value
cvodes.CVodeGetNumNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumNonlinSolvConvFails.restype = ctypes.c_int

def CVodeGetNonlinSolvStats(cvodememobj):
	"""A convenience function that provides the nonlinear solver optional outputs in a tuple (NumNonlinSolvIters, NumNonlinsolvConvFails).
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nniters = ctypes.c_long()
	nncfails = ctypes.c_long()
	ret = cvodes.CVodeGetNonlinSolvStats(cvodememobj.obj, ctypes.byref(nniters), ctypes.byref(nncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNonlinSolvStats() failed with flag %i"%(ret))
	return (nniters.value, nncfails.value)
cvodes.CVodeGetNonlinSolvStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNonlinSolvStats.restype = ctypes.c_int

def CVodeGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVODE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVodeGetReturnFlagName(flag)
cvodes.CVodeGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVodeGetReturnFlagName.restype = ctypes.c_char_p

def CVodeGetQuad(cvodememobj, t, yQout):
	ret = cvodes.CVodeGetQuad(cvodememobj.obj, t, yQout.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuad() failed with flag %i"%(ret))
cvodes.CVodeGetQuad.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetQuad.restype = ctypes.c_int

def CVodeGetQuadDky(cvodememobj, t, k, dky):
	ret = cvodes.CVodeGetQuadDky(cvodememobj.obj, t, k, dky.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadDky() failed with flag %i"%(ret))
cvodes.CVodeGetQuadDky.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetQuadDky.restype = ctypes.c_int

def CVodeGetQuadNumRhsEvals(cvodememobj):
	nfQevals = ctypes.c_long(0)
	ret = cvodes.CVodeGetQuadNumRhsEvals(cvodememobj.obj, ctypes.byref(nfQevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadNumRhsEvals() failed with flag %i"%(ret))
	return nfQevals.value
cvodes.CVodeGetQuadNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetQuadNumRhsEvals.restype = ctypes.c_int

def CVodeGetQuadNumErrTestFails(cvodememobj):
	nQetfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetQuadNumErrTestFails(cvodememobj.obj, ctypes.byref(nQetfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadNumErrTestFails() failed with flag %i"%(ret))
	return nQetfails.value
cvodes.CVodeGetQuadNumErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetQuadNumErrTestFails.restype = ctypes.c_int

def CVodeGetQuadErrWeights(cvodememobj, eQweight):
	ret = cvodes.CVodeGetQuadErrWeights(cvodememobj.obj, eQweight.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadErrWeights() failed with flag %i"%(ret))
cvodes.CVodeGetQuadErrWeights.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetQuadErrWeights.restype = ctypes.c_int

def CVodeGetQuadStats(cvodememobj):
	nfQevals = ctypes.c_long(0)
	nQetfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetQuadStats(cvodememobj.obj, ctypes.byref(nfQevals), ctypes.byref(nQetfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadStats() failed with flag %i"%(ret))
	return (nfQevals.value, nQetfails.value)
cvodes.CVodeGetQuadStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetQuadStats.restype = ctypes.c_int

def CVodeGetSens(cvodememobj, t, ySout):
	ret = cvodes.CVodeGetSens(cvodememobj.obj, t, ySout.cdata)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSens() failed with flag %i"%(ret))
cvodes.CVodeGetSens.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(ctypes.POINTER(nvecserial._NVector))]
cvodes.CVodeGetSens.restype = ctypes.c_int

def CVodeGetSens1(cvodememobj, t, i, ySout):
	ret = cvodes.CVodeGetSens1(cvodememobj.obj, t, i, ySout.cdata)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSens1() failed with flag %i"%(ret))
cvodes.CVodeGetSens1.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetSens1.restype = ctypes.c_int

def CVodeGetSensDky(cvodememobj):
	t = ctypes.c_int(0)
	ret = cvodes.CVodeGetSensDky(cvodememobj.obj, ctypes.byref(t), k, dkyA.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSensDky() failed with flag %i"%(ret))
	return t.value
cvodes.CVodeGetSensDky.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(nvecserial._NVector))]
cvodes.CVodeGetSensDky.restype = ctypes.c_int

def CVodeGetSensDky1(cvodememobj, t, k, i, dky):
	ret = cvodes.CVodeGetSensDky1(cvodememobj.obj, t, k, i, dky.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSensDky1() failed with flag %i"%(ret))
cvodes.CVodeGetSensDky1.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetSensDky1.restype = ctypes.c_int

def CVodeGetNumSensRhsEvals(cvodememobj):
	nfSevals = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumSensRhsEvals(cvodememobj.obj, ctypes.byref(nfSevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSensRhsEvals() failed with flag %i"%(ret))
	return nfSevals.value
cvodes.CVodeGetNumSensRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSensRhsEvals.restype = ctypes.c_int

def CVodeGetNumRhsEvalsSens(cvodememobj):
	nfevalsS = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumRhsEvalsSens(cvodememobj.obj, ctypes.byref(nfevalsS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumRhsEvalsSens() failed with flag %i"%(ret))
	return nfevalsS.value
cvodes.CVodeGetNumRhsEvalsSens.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumRhsEvalsSens.restype = ctypes.c_int

def CVodeGetNumSensErrTestFails(cvodememobj):
	nSetfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumSensErrTestFails(cvodememobj.obj, ctypes.byref(nSetfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSensErrTestFails() failed with flag %i"%(ret))
	return nSetfails.value
cvodes.CVodeGetNumSensErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSensErrTestFails.restype = ctypes.c_int

def CVodeGetNumSensLinSolvSetups(cvodememobj):
	nlinsetupsS = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumSensLinSolvSetups(cvodememobj.obj, ctypes.byref(nlinsetupsS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSensLinSolvSetups() failed with flag %i"%(ret))
	return nlinsetupsS.value
cvodes.CVodeGetNumSensLinSolvSetups.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSensLinSolvSetups.restype = ctypes.c_int

def CVodeGetSensErrWeights(cvodememobj, eSweight):
	ret = cvodes.CVodeGetSensErrWeights(cvodememobj.obj, eSweight.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSensErrWeights() failed with flag %i"%(ret))
cvodes.CVodeGetSensErrWeights.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetSensErrWeights.restype = ctypes.c_int

def CVodeGetSensStats(cvodememobj):
	nfSevals = ctypes.c_long(0)
	nfevalsS = ctypes.c_long(0)
	nSetfails = ctypes.c_long(0)
	nlinsetupS = ctypes.c_long(0)
	ret = cvodes.CVodeGetSensStats(cvodememobj.obj, ctypes.byref(nfSevals), ctypes.byref(nfevalsS), ctypes.byref(nSetfails), ctypes.byref(nlinsetupsS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSensStats() failed with flag %i"%(ret))
	return (nfSevals.value, nfevalsS.value, nSetfails.value, nlinsetupS.value)
cvodes.CVodeGetSensStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetSensStats.restype = ctypes.c_int

def CVodeGetNumSensNonlinSolvIters(cvodememobj):
	nSniters = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumSensNonlinSolvIters(cvodememobj.obj, ctypes.byref(nSniters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSensNonlinSolvIters() failed with flag %i"%(ret))
	return nSniters.value
cvodes.CVodeGetNumSensNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSensNonlinSolvIters.restype = ctypes.c_int

def CVodeGetNumSensNonlinSolvConvFails(cvodememobj):
	nSncfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumSensNonlinSolvConvFails(cvodememobj.obj, ctypes.byref(nSncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSensNonlinSolvConvFails() failed with flag %i"%(ret))
	return nSncfails.value
cvodes.CVodeGetNumSensNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumSensNonlinSolvConvFails.restype = ctypes.c_int

def CVodeGetNumStgrSensNonlinSolvIters(cvodememobj):
	nSTGR1niters = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumStgrSensNonlinSolvIters(cvodememobj.obj, ctypes.byref(nSTGR1niters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumStgrSensNonlinSolvIters() failed with flag %i"%(ret))
	return nSTGR1niters.value
cvodes.CVodeGetNumStgrSensNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumStgrSensNonlinSolvIters.restype = ctypes.c_int

def CVodeGetNumStgrSensNonlinSolvConvFails(cvodememobj):
	nSTGR1ncfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetNumStgrSensNonlinSolvConvFails(cvodememobj.obj, ctypes.byref(nSTGR1ncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumStgrSensNonlinSolvConvFails() failed with flag %i"%(ret))
	return nSTGR1ncfails.value
cvodes.CVodeGetNumStgrSensNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetNumStgrSensNonlinSolvConvFails.restype = ctypes.c_int

def CVodeGetSensNonlinSolvStats(cvodememobj):
	nSniters = ctypes.c_long(0)
	nSncfails = ctypes.c_long(0)
	ret = cvodes.CVodeGetSensNonlinSolvStats(cvodememobj.obj, ctypes.byref(nSniters), ctypes.byref(nSncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetSensNonlinSolvStats() failed with flag %i"%(ret))
	return (nSniters.value, nSncfails.value)
cvodes.CVodeGetSensNonlinSolvStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVodeGetSensNonlinSolvStats.restype = ctypes.c_int

#def CVodeFree(cvodememobj):
#	"""CVodeFree frees the problem memory cvodememobj allocated by CVodeCreate and CVodeMalloc."""
#	cvodes.CVodeFree(cvodememobj.obj)
#	del cvodememobj.obj

def CVodeQuadFree(cvodememobj):
	ret = cvodes.CVodeQuadFree(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeQuadFree() failed with flag %i"%(ret))
cvodes.CVodeQuadFree.argtypes = [ctypes.c_void_p]
cvodes.CVodeQuadFree.restype = None

def CVodeSensFree(cvodememobj):
	ret = cvodes.CVodeSensFree(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSensFree() failed with flag %i"%(ret))
cvodes.CVodeSensFree.argtypes = [ctypes.c_void_p]
cvodes.CVodeSensFree.restype = None

def CVadjMalloc(cvodememobj, steps, interp):
	ret = cvodes.CVadjMalloc(cvodememobj.obj, steps, interp)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjMalloc() failed with flag %i"%(ret))
	return ret
cvodes.CVadjMalloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_int]
cvodes.CVadjMalloc.restype = ctypes.c_void_p

def CVadjSetInterpType(cvadj_mem, interp):
	ret = cvodes.CVadjSetInterpType(cvadj_mem, interp)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjSetInterpType() failed with flag %i"%(ret))
cvodes.CVadjSetInterpType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVadjSetInterpType.restype = ctypes.c_int

def CVodeF(cvadj_mem, tout, yout, tret, itask, ncheckPtr):
	ret = cvodes.CVodeF(cvadj_mem, tout, yout.data, tret, itask, ncheckPtr)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeF() failed with flag %i"%(ret))
cvodes.CVodeF.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype ), ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
cvodes.CVodeF.restype = ctypes.c_int

def CVodeCreateB(cvadj_mem, lmmB, iterB):
	ret = cvodes.CVodeCreateB(cvadj_mem, lmmB, iterB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeCreateB() failed with flag %i"%(ret))
cvodes.CVodeCreateB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVodeCreateB.restype = ctypes.c_int

def CVodeMallocB(cvadj_mem, fB, tB0, yB0, itolB, reltolB, abstolB):
	if itolB == CV_SS:
		if type(abstolB) == realtype:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, ctypes.byref(abstolB))
		elif type(abstolB) == float or type(abstolB) == int:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, ctypes.byref(realtype(abstolB)))
		elif abstolB == None:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be a floating point number if itolB is CV_SS")
	elif itolB == CV_SV:
		if type(abstolB) == NVector:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB.data)
		elif abstolB == None:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be an NVector if itolB is CV_SV")
	elif itolB == CV_WF:
		if abstolB == None:
			ret = cvodes.CVodeMallocB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be None if itolB is CV_WF")
	else:
		raise ValueError("itolB must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeMallocB() failed with flag %i"%(ret))
cvodes.CVodeMallocB.argtypes = [ctypes.c_void_p, CVRhsFnB, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeMallocB.restype = ctypes.c_int

def CVodeSetErrHandlerFnB(cvadj_mem, ehfunB, eh_dataB):
	ret = cvodes.CVodeSetErrHandlerFnB(cvadj_mem, WrapCallbackCVErrHandlerFn(ehfunB), eh_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrHandlerFnB() failed with flag %i"%(ret))
cvodes.CVodeSetErrHandlerFnB.argtypes = [ctypes.c_void_p, CVErrHandlerFn, ctypes.c_void_p]
cvodes.CVodeSetErrHandlerFnB.restype = ctypes.c_int

def CVodeSetErrFileB(cvadj_mem, fileobj):
	ret = cvodes.CVodeSetErrFileB(cvadj_mem, sundials_core.fdopen(fileobj.fileno, filobj.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrFileB() failed with flag %i"%(ret))
cvodes.CVodeSetErrFileB.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetErrFileB.restype = ctypes.c_int

def CVodeSetIterTypeB(cvadj_mem, iterB):
	ret = cvodes.CVodeSetIterTypeB(cvadj_mem, iterB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetIterTypeB() failed with flag %i"%(ret))
cvodes.CVodeSetIterTypeB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetIterTypeB.restype = ctypes.c_int

def CVodeSetFdataB(cvadj_mem, f_dataB):
	ret = cvodes.CVodeSetFdataB(cvadj_mem, f_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetFdataB() failed with flag %i"%(ret))
cvodes.CVodeSetFdataB.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetFdataB.restype = ctypes.c_int

def CVodeSetMaxOrdB(cvadj_mem, maxordB):
	ret = cvodes.CVodeSetMaxOrdB(cvadj_mem, maxordB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxOrdB() failed with flag %i"%(ret))
cvodes.CVodeSetMaxOrdB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetMaxOrdB.restype = ctypes.c_int

def CVodeSetMaxNumStepsB(cvadj_mem, mxstepsB):
	ret = cvodes.CVodeSetMaxNumStepsB(cvadj_mem, mxstepsB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNumStepsB() failed with flag %i"%(ret))
cvodes.CVodeSetMaxNumStepsB.argtypes = [ctypes.c_void_p, ctypes.c_long]
cvodes.CVodeSetMaxNumStepsB.restype = ctypes.c_int

def CVodeSetStabLimDetB(cvadj_mem, stldetB):
	ret = cvodes.CVodeSetStabLimDetB(cvadj_mem, stldetB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStabLimDetB() failed with flag %i"%(ret))
cvodes.CVodeSetStabLimDetB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVodeSetStabLimDetB.restype = ctypes.c_int

def CVodeSetInitStepB(cvadj_mem, hinB):
	ret = cvodes.CVodeSetInitStepB(cvadj_mem, hinB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetInitStepB() failed with flag %i"%(ret))
cvodes.CVodeSetInitStepB.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetInitStepB.restype = ctypes.c_int

def CVodeSetMinStepB(cvadj_mem, hminB):
	ret = cvodes.CVodeSetMinStepB(cvadj_mem, hminB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMinStepB() failed with flag %i"%(ret))
cvodes.CVodeSetMinStepB.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetMinStepB.restype = ctypes.c_int

def CVodeSetMaxStepB(cvadj_mem, hmaxB):
	ret = cvodes.CVodeSetMaxStepB(cvadj_mem, hmaxB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxStepB() failed with flag %i"%(ret))
cvodes.CVodeSetMaxStepB.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVodeSetMaxStepB.restype = ctypes.c_int

def CVodeReInitB(cvadj_mem, fB, tB0, yB0, itolB, reltolB, abstolB):
	if itolB == CV_SS:
		if type(abstolB) == realtype:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, ctypes.byref(abstolB))
		elif type(abstolB) == float or type(abstolB) == int:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, ctypes.byref(realtype(abstolB)))
		elif abstolB == None:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be a floating point number if itolB is CV_SS")
	elif itolB == CV_SV:
		if type(abstolB) == NVector:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB.data)
		elif abstolB == None:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be an NVector if itolB is CV_SV")
	elif itolB == CV_WF:
		if abstolB == None:
			ret = cvodes.CVodeReInitB(cvadj_mem, WrapCallbackCVRhsFnB(fB), tB0, yB0.data, itolB, reltolB, abstolB)
		else:
			raise TypeError("abstolB must be None if itolB is CV_WF")
	else:
		raise ValueError("itolB must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeReInitB() failed with flag %i"%(ret))
cvodes.CVodeReInitB.argtypes = [ctypes.c_void_p, CVRhsFnB, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeReInitB.restype = ctypes.c_int

def CVodeSetQuadFdataB(cvadj_mem, fQ_dataB):
	ret = cvodes.CVodeSetQuadFdataB(cvadj_mem, fQ_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetQuadFdataB() failed with flag %i"%(ret))
cvodes.CVodeSetQuadFdataB.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvodes.CVodeSetQuadFdataB.restype = ctypes.c_int

def CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB):
	if itolQB == CV_SS:
		if type(abstolQB) == realtype:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, ctypes.byref(abstolQB))
		elif type(abstolQB) == float or type(abstolQB) == int:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, ctypes.byref(realtype(abstolQB)))
		elif abstolQB == None:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB)
		else:
			raise TypeError("abstolQB must be a floating point number if itolQB is CV_SS")
	elif itolQB == CV_SV:
		if type(abstolQB) == NVector:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB.data)
		elif abstolQB == None:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB)
		else:
			raise TypeError("abstolQB must be an NVector if itol is CV_SV")
	elif itolQB == CV_WF:
		if abstolQB == None:
			ret = cvodes.CVodeSetQuadErrConB(cvadj_mem, errconQB, itolQB, reltolQB, abstolQB)
		else:
			raise TypeError("abstolQB must be None if itolQB is CV_WF")
	else:
		raise ValueError("itolQB must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetQuadErrConB() failed with flag %i"%(ret))
cvodes.CVodeSetQuadErrConB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, realtype, ctypes.c_void_p]
cvodes.CVodeSetQuadErrConB.restype = ctypes.c_int

def CVodeQuadMallocB(cvadj_mem, fQB, yQB0):
	ret = cvodes.CVodeQuadMallocB(cvadj_mem, WrapCallbackCVQuadRhsFnB(fQB), yQB0.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeQuadMallocB() failed with flag %i"%(ret))
cvodes.CVodeQuadMallocB.argtypes = [ctypes.c_void_p, CVQuadRhsFnB, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeQuadMallocB.restype = ctypes.c_int

def CVodeQuadReInitB(cvadj_mem, fQB, yQB0):
	ret = cvodes.CVodeQuadReInitB(cvadj_mem, WrapCallbackCVQuadRhsFnB(fQB), yQB0.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeQuadReInitB() failed with flag %i"%(ret))
cvodes.CVodeQuadReInitB.argtypes = [ctypes.c_void_p, CVQuadRhsFnB, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeQuadReInitB.restype = ctypes.c_int

def CVodeB(cvadj_mem, tBout, yBout, tBret, itaskB):
	ret = cvodes.CVodeB(cvadj_mem, tBout, yBout.data, tBret, itaskB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeB() failed with flag %i"%(ret))
cvodes.CVodeB.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype ), ctypes.c_int]
cvodes.CVodeB.restype = ctypes.c_int

def CVodeGetQuadB(cvadj_mem, qB):
	ret = cvodes.CVodeGetQuadB(cvadj_mem, qB.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetQuadB() failed with flag %i"%(ret))
cvodes.CVodeGetQuadB.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVodeGetQuadB.restype = ctypes.c_int

def CVadjFree(cvadj_mem):
	ret = cvodes.CVadjFree(cvadj_mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjFree() failed with flag %i"%(ret))
cvodes.CVadjFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvodes.CVadjFree.restype = None

def CVadjGetCVodeBmem(cvadj_mem):
	ret = cvodes.CVadjGetCVodeBmem(cvadj_mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetCVodeBmem() failed with flag %i"%(ret))
	return CVodeMemObj(ret) 
cvodes.CVadjGetCVodeBmem.argtypes = [ctypes.c_void_p]
cvodes.CVadjGetCVodeBmem.restype = ctypes.c_void_p

def CVadjGetReturnFlagName(flag):
	return cvodes.CVadjGetReturnFlagName(flag)
cvodes.CVadjGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVadjGetReturnFlagName.restype = ctypes.c_char_p

def CVadjGetY(cvadj_mem, t, y):
	ret = cvodes.CVadjGetY(cvadj_mem, t, y.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetY() failed with flag %i"%(ret))
cvodes.CVadjGetY.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector)]
cvodes.CVadjGetY.restype = ctypes.c_int

class CVadjCheckPointRec(ctypes.Structure):
	_fields_ = [
		("my_addr", ctypes.c_void_p),
		("next_addr", ctypes.c_void_p),
		("t0", realtype),
		("t1", realtype),
		("nstep", ctypes.c_long),
		("order", ctypes.c_int),
		("step", realtype)
	]

def CVadjGetCheckPointsInfo(cvadj_mem):
	ckpnt = CVadjCheckPointRec()
	ret = cvodes.CVadjGetCheckPointsInfo(cvadj_mem, ctypes.byref(ckpnt))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetCheckPointsInfo() failed with flag %i"%(ret))
	return ckpnt
cvodes.CVadjGetCheckPointsInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(CVadjCheckPointRec)]
cvodes.CVadjGetCheckPointsInfo.restype = ctypes.c_int

def CVadjGetDataPointHermite(cvadj_memi, which, t, y, yd):
	ret = cvodes.CVadjGetDataPointHermite(cvadj_mem, which, ctypes.byref(t), y.data, yd.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetDataPointHermite() failed with flag %i"%(ret))
cvodes.CVadjGetDataPointHermite.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(realtype ), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector)]
cvodes.CVadjGetDataPointHermite.restype = ctypes.c_int

def CVadjGetDataPointPolynomial(cvadj_mem, which, t, order, y):
	ret = cvodes.CVadjGetDataPointPolynomial(cvadj_mem, which, ctypes.byref(t), ctypes.byref(order), y.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetDataPointPolynomial() failed with flag %i"%(ret))
cvodes.CVadjGetDataPointPolynomial.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(realtype ), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(nvecserial._NVector)]
cvodes.CVadjGetDataPointPolynomial.restype = ctypes.c_int

def CVadjGetCurrentCheckPoint(cvadj_mem):
	addr = ctypes.c_int(0)
	ret = cvodes.CVadjGetCurrentCheckPoint(cvadj_mem, ctypes.byref(addr))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVadjGetCurrentCheckPoint() failed with flag %i"%(ret))
	return addr.value
cvodes.CVadjGetCurrentCheckPoint.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_void_p)]
cvodes.CVadjGetCurrentCheckPoint.restype = ctypes.c_int

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

The callable python object must take exactly 3 parameters, which will be passed in as
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
	ret = cvodes.ModifiedGS(v.data, h, k, p, new_vk_norm)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ModifiedGS() failed with flag %i"%(ret))
cvodes.ModifiedGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype)]
cvodes.ModifiedGS.restype = ctypes.c_int

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
	ret = cvodes.ClassicalGS(v.data, h, k, p, new_vk_norm, temp.data, s)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ClassicalGS() failed with flag %i"%(ret))
cvodes.ClassicalGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype)]
cvodes.ClassicalGS.restype = ctypes.c_int

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
	ret = cvodes.QRfact(n, h, q, job)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRfact() failed with flag %i"%(ret))
cvodes.QRfact.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.c_int]
cvodes.QRfact.restype = ctypes.c_int

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
	ret = cvodes.QRsol(n, h, q, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRsol() failed with flag %i"%(ret))
cvodes.QRsol.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
cvodes.QRsol.restype = ctypes.c_int


#########################
# sundials_smalldense.h #
#########################

def denalloc(m, n):
	"""Allocates storage for an m by n small dense matrix and returns a pointer to the newly allocated storage if successful. If the memory request cannot be satisfied, a ValueError Exception is raised.

The underlying type of the dense matrix returned is a double pointer to realtype. If we allocate a dense matrix a by a = denalloc(m, n), then a[j][i] references the (i,j)-th element of the matrix a, 0 <= i < m, 0 <= j < n,  and a[j] is a pointer to the first element in the jth column of a. The location a[0] contains a pointer to m*n contiguous locations which contain the elements of a."""
	ret = cvodes.denalloc(m, n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denalloc could not allocate memory')
cvodes.denalloc.argtypes = [ctypes.c_long, ctypes.c_long]
cvodes.denalloc.restype = ctypes.POINTER(ctypes.POINTER(realtype))

def denallocpiv(n):
	"""Allocates an array of n long int, and returns a pointer to the first element in the array if successful, otherwise a ValueError exception is raised."""
	ret = cvodes.denallocpiv(n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denallocpiv could not allocate memory')
cvodes.denallocpiv.argtypes = [ctypes.c_long]
cvodes.denallocpiv.restype = ctypes.POINTER(ctypes.c_long)

def denGETRF(a, m, n, p):
	"""Factors the m by n dense matrix a, m>=n. It overwrites the elements of a with its LU factors and keeps track of the pivot rows chosen in the pivot array p. A successful LU factorization leaves the matrix a and the pivot array p with the following information:
	(1) p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., n-1.
	(2) If the unique LU factorization of a is given by Pa = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1.0 on the diagonal, and U is an upper triangular matrix, then the upper triangular part of a (including its diagonal) contains U and the strictly lower trapezoidal part of a contains the multipliers, I-L.

Note that for square matrices (m=n), L is unit lower triangular.

denGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	return cvodes.denGETRF(a, m, n, p)
cvodes.denGETRF.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]
cvodes.denGETRF.restype = ctypes.c_long

def denGETRS(a, n, p, b):
	"""Solves the n by n linear system a*x = b. It assumes that a has been LU factored and the pivot array p has been set by a successful call to denGETRF(a,n,n,p).
	
denGETRS does not check whether a is square! The solution x is written into the b array."""
	return cvodes.denGETRS(a, n, p, b)
cvodes.denGETRS.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
cvodes.denGETRS.restype = None

def denzero(a, m, n):
	"""Sets all the elements of the m by n dense matrix a to be 0.0."""
	cvodes.denzero(a, m, n)
cvodes.denzero.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvodes.denzero.restype = None

def dencopy(a, b, m, n):
	"""Copies the m by n dense matrix a into the m by n dense matrix b."""
	cvodes.dencopy(a, b, m, n)
cvodes.dencopy.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvodes.dencopy.restype = None

def denscale(c, a, m, n):
	"""Scales every element in the m by n dense matrix a by c."""
	cvodes.denscale(c, a, m, n)
cvodes.denscale.argtypes = [realtype, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvodes.denscale.restype = None

def denaddI(a, n):
	"""Increments the diagonal elements of the dense m by n matrix a by 1.0. (a_ii <= a_ii + 1, i=1,2,..n-1.)
	
denaddI is typically used with square matrices.
denaddI does NOT check for m >= n! Therefore, a segmentation fault will occur if m<n!
"""
	cvodes.denaddI(a, n)
cvodes.denaddI.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long]
cvodes.denaddI.restype = None

def denfreepiv(p):
	"""Frees the pivot array p allocated by denallocpiv."""
	cvodes.denfreepiv(p)
cvodes.denfreepiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvodes.denfreepiv.restype = None

def denfree(a):
	"""Frees the dense matrix a allocated by denalloc."""
	cvodes.denfree(a)
cvodes.denfree.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype))]
cvodes.denfree.restype = None

def denprint(a, m, n):
	"""Prints the m by n dense matrix a to standard output as it would normally appear on paper. It is intended as a debugging tool with small values of m and n. The elements are printed using the %g/%lg/%Lg option. A blank line is printed before and after the matrix."""
	cvodes.denprint(a, m, n)
cvodes.denprint.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvodes.denprint.restype = None

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
	return cvodes.SpbcgMalloc(l_max, vec_tmpl.data)
cvodes.SpbcgMalloc.argtypes = [ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvodes.SpbcgMalloc.restype = SpbcgMem

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
	ret = cvodes.SpbcgSolve(mem, A_data, x.data, b.data, pretype, delta, P_data, sx.data, sb.data, WrapCallbackATimesFn(atimes), WrapCallbackPSolveFn(psolve), res_norm, nli, nps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgSolve() failed with flag %i"%(ret))
cvodes.SpbcgSolve.argtypes = [SpbcgMem, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ATimesFn, PSolveFn, ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
cvodes.SpbcgSolve.restype = ctypes.c_int

def SpbcgFree(mem):
	"""Deallocates any memory allocated implicitly by Spbcg"""
	ret = cvodes.SpbcgFree(mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgFree() failed with flag %i"%(ret))
cvodes.SpbcgFree.argtypes = [SpbcgMem]
cvodes.SpbcgFree.restype = None


#################
# cvodes_diag.h #
#################

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

#Return values for adjoint module
CVDIAG_ADJMEM_NULL = -101

def CVDiag(cvodememobj):
	"""Links the main integrator with the CVDIAG linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	ret = cvodes.CVDiag(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiag() failed with flag %i"%(ret))
cvodes.CVDiag.argtypes = [ctypes.c_void_p]
cvodes.CVDiag.restype = ctypes.c_int

def CVDiagGetWorkSpace(cvodememobj):
	"""Returns the CVDIAG real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvodes.CVDiagGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVDiagGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVDiagGetWorkSpace.restype = ctypes.c_int

def CVDiagGetNumRhsEvals(cvodememobj):
	"""Returns the number of calls to the user f routine due to finite difference Jacobian evaluation.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()

Note: The number of diagonal approximate Jacobians formed is equal to the number of CVDiagSetup calls. This number is available through CVodeGetNumLinSolvSetups."""
	retval = ctypes.c_long(0)
	ret = cvodes.CVDiagGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVDiagGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVDiagGetNumRhsEvals.restype = ctypes.c_int

def CVDiagGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVDIAG interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int(0)
	ret = cvodes.CVDiagGetLastFlag(cvodes_mem, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetLastFlag() failed with flag %i"%(ret))
	return retval
cvodes.CVDiagGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVDiagGetLastFlag.restype = ctypes.c_int

def CVDiagGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDIAG return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVDiagGetReturnFlagName(flag)
cvodes.CVDiagGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVDiagGetReturnFlagName.restype = ctypes.c_char_p

def CVDiagB(cvadj_mem):
	"""CVDiagB links the main CVODES integrator with the CVDIAG linear solver for the backward integration."""
	ret = cvodes.CVDiagB(cvadj_mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagB() failed with flag %i"%(ret))
cvodes.CVDiagB.argtypes = [ctypes.c_void_p]
cvodes.CVDiagB.restype = ctypes.c_int


###################
# cvodes_bbdpre.h #
###################

CVBBDPRE_SUCCESS = 0
CVBBDPRE_PDATA_NULL = -11
CVBBDPRE_FUNC_UNRECVR = -12

CVBBDPRE_ADJMEM_NULL = -111
CVBBDPRE_PDATAB_NULL = -112
CVBBDPRE_MEM_FAIL = -113

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
	ret = cvodes.CVBBDPrecAlloc(cvodememobj.obj, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, WrapCallbackCVLocalFn(gloc), WrapCallbackCVCommFn(cfn))
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecAlloc() failed to allocate memory")
	return ret
cvodes.CVBBDPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, CVLocalFn, CVCommFn]
cvodes.CVBBDPrecAlloc.restype = ctypes.c_void_p

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
	ret = cvodes.CVBBDSptfqmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSptfqmr() failed with flag %i"%(ret))
cvodes.CVBBDSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBBDSptfqmr.restype = ctypes.c_int

def CVBBDSpbcg(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSptfqmr links the CVBBDPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions:
	1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory;
	2) Sets the preconditioner data structure for CVSPTFQMR
	3) Sets the preconditioner setup routine for CVSPTFQMR
	4) Sets the preconditioner solve routine for CVSPTFQMR

Its first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSptfqmr."""
	ret = cvodes.CVBBDSpbcg(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpbcg() failed with flag %i"%(ret))
cvodes.CVBBDSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBBDSpbcg.restype = ctypes.c_int

def CVBBDSpgmr(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSpbcg links the CVBBDPRE preconditioner to the CVSPBCG linear solver. It performs the following actions:
	1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory;
	2) Sets the preconditioner data structure for CVSPBCG
	3) Sets the preconditioner setup routine for CVSPBCG
	4) Sets the preconditioner solve routine for CVSPBCG

Its first 3 arguments are the same as for CVSpbcg (see cvspbcg.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSpbcg."""
	ret = cvodes.CVBBDSpgmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpgmr() failed with flag %i"%(ret))
cvodes.CVBBDSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBBDSpgmr.restype = ctypes.c_int

def CVBBDPrecReInit(bbd_data, mudq, mldq, dqrely, gloc, cfn):
	"""CVBBDPrecReInit re-initializes the BBDPRE module when solving a sequence of problems of the same size with CVSPGMR/CVBBDPRE or CVSPBCG/CVBBDPRE or CVSPTFQMR/CVBBDPRE provided there is no change in Nlocal, mukeep, or mlkeep. After solving one problem, and after calling CVodeReInit to re-initialize the integrator for a subsequent problem, call CVBBDPrecReInit. Then call CVSpgmrSet* or CVSpbcgSet* or CVSptfqmrSet* functions if necessary for any changes to CVSPGMR, CVSPBCG, or CVSPTFQMR parameters, before calling CVode.

The first argument to CVBBDPrecReInit must be the pointer pdata that was returned by CVBBDPrecAlloc. All other arguments have the same names and meanings as those of CVBBDPrecAlloc."""
	ret = cvodes.CVBBDPrecReInit(bbd_data, mudq, mldq, dqrely, WrapCallbackCVLocalFn(gloc), WrapCallbackCVCommFn(cfn))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecReInit() failed with flag %i"%(ret))
cvodes.CVBBDPrecReInit.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, realtype, CVLocalFn, CVCommFn]
cvodes.CVBBDPrecReInit.restype = ctypes.c_int

def CVBBDPrecFree(bbd_data):
	"""CVBBDPrecFree frees the memory block bbd_data allocated by the call to CVBBDAlloc."""
	ret = cvodes.CVBBDPrecFree(bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecFree() failed with flag %i"%(ret))
cvodes.CVBBDPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvodes.CVBBDPrecFree.restype = None

def CVBBDPrecGetWorkSpace(bbd_data):
	"""CVBBDPrecGetWorkSpace returns the CVBBDPrec real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvodes.CVBBDPrecGetWorkSpace(bbd_data, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVBBDPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVBBDPrecGetWorkSpace.restype = ctypes.c_int

def CVBBDPrecGetNumGfnEvals(bbd_data):
 	"""CVBBDPrecGetNumGfnEvals returns the number of calls to gfn."""
	ngevalsBBDP = ctypes.c_long(0)
	ret = cvodes.CVBBDPrecGetNumGfnEvals(bbd_data, ctypes.byref(ngevalsBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecGetNumGfnEvals() failed with flag %i"%(ret))
	return ngevalsBBDP.value
cvodes.CVBBDPrecGetNumGfnEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVBBDPrecGetNumGfnEvals.restype = ctypes.c_int

def CVBBDPrecGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBBDPRE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVBBDPrecGetReturnFlagName(flag)
cvodes.CVBBDPrecGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVBBDPrecGetReturnFlagName.restype = ctypes.c_char_p

CVLocalFnB = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVLocalFnB(func):
	exec 'def __CallbackInterface_%s(NlocalB, t, y, yB, gB, f_dataB):\n\treturn __ActualCallback[%i](NlocalB, t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(gB), f_dataB)'%(func.func_name, len(__ActualCallback))
	if func == None:
		return ctypes.cast(None, CVLocalFnB)
	__ActualCallback.append(func)
	tmp = CVLocalFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVCommFnB = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVCommFnB(func):
	exec 'def __CallbackInterface_%s(NlocalB, t, y, yB, f_dataB):\n\treturn __ActualCallback[%i](NlocalB, t, nvecserial.NVector(y), nvecserial.NVector(yB), f_dataB)'%(func.func_name, len(__ActualCallback))
	if func == None:
		return ctypes.cast(None, CVCommFnB)
	__ActualCallback.append(func)
	tmp = CVCommFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVBBDPrecAllocB(cvadj_mem, NlocalB, mudqB, mldqB, mukeepB, mlkeepB, dqrelyB, glocB, cfnB):
	ret = cvodes.CVBBDPrecAllocB(cvadj_mem, NlocalB, mudqB, mldqB, mukeepB, mlkeepB, dqrelyB, WrapCallbackCVLocalFnB(glocB), WrapCallbackCVCommFnB(cfnB))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecAllocB() failed with flag %i"%(ret))
cvodes.CVBBDPrecAllocB.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, CVLocalFnB, CVCommFnB]
cvodes.CVBBDPrecAllocB.restype = ctypes.c_int

def CVBBDSptfqmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBBDSptfqmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSptfqmrB() failed with flag %i"%(ret))
cvodes.CVBBDSptfqmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBBDSptfqmrB.restype = ctypes.c_int

def CVBBDSpbcgB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBBDSpbcgB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpbcgB() failed with flag %i"%(ret))
cvodes.CVBBDSpbcgB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBBDSpbcgB.restype = ctypes.c_int

def CVBBDSpgmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBBDSpgmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpgmrB() failed with flag %i"%(ret))
cvodes.CVBBDSpgmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBBDSpgmrB.restype = ctypes.c_int

def CVBBDPrecReInitB(cvadj_mem, mudqB, mldqB, dqrelyB, glocB, cfnB):
	ret = cvodes.CVBBDPrecReInitB(cvadj_mem, mudqB, mldqB, dqrelyB, WrapCallbackCVLocalFnB(glocB), WrapCallbackCVCommFnB(cfnB))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecReInitB() failed with flag %i"%(ret))
cvodes.CVBBDPrecReInitB.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, realtype, CVLocalFnB, CVCommFnB]
cvodes.CVBBDPrecReInitB.restype = ctypes.c_int

def CVBBDPrecFreeB(cvadj_mem):
	ret = cvodes.CVBBDPrecFreeB(cvadj_mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecFreeB() failed with flag %i"%(ret))
cvodes.CVBBDPrecFreeB.argtypes = [ctypes.c_void_p]
cvodes.CVBBDPrecFreeB.restype = None

####################
# cvodes_bandpre.h #
####################

CVBANDPRE_SUCCESS = 0
CVBANDPRE_PDATA_NULL = -11
CVBANDPRE_RHSFUNC_UNRECVR = -12

CVBANDPRE_ADJMEM_NULL = -111
CVBANDPRE_MEM_FAIL = -112

def CVBandPrecAlloc(cvodememobj, N, mu, ml):
	"""CVBandPrecAlloc allocates and initializes a CVBandPrecData structure to be passed to CVSp* (and subsequently used by CVBandPrecSetup and CVBandPrecSolve). The parameters of CVBandPrecAlloc are as follows:
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate() 
	N	the problem size.  
	mu	the upper half bandwidth.  
	ml	the lower half bandwidth.
CVBandPrecAlloc returns the storage pointer of type CVBandPrecData, or NULL if the request for storage cannot be satisfied.  
NOTE:	The band preconditioner assumes a serial implementation of the NVECTOR package. Therefore, CVBandPrecAlloc will first test for a compatible N_Vector internal representation by checking for required functions."""
	ret = cvodes.CVBandPrecAlloc(cvodememobj.obj, N, mu, ml)
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecAlloc() failed to allocate memory")
	return ret
cvodes.CVBandPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvodes.CVBandPrecAlloc.restype = ctypes.c_void_p

def CVBPSptfqmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSptfqmr links the CVBANDPPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions: 
	1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPTFQMR 
	3) Sets the preconditioner setup routine for CVSPTFQMR 
	4) Sets the preconditioner solve routine for CVSPTFQMR 

	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done. This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPTFQMR solver. Pass 0 to use the default value CVSPILS_MAXL=5."""
	ret = cvodes.CVBPSptfqmr(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSptfqmr() failed with flag %i"%(ret))
cvodes.CVBPSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBPSptfqmr.restype = ctypes.c_int

def CVBPSpbcg(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpbcg links the CVBANDPPRE preconditioner to the CVSPBCG linear solver. It performs the following actions: 
	1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPBCG 
	3) Sets the preconditioner setup routine for CVSPBCG 
	4) Sets the preconditioner solve routine for CVSPBCG 

	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in iterative.h. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPBCG solver. Pass 0 to use the default value CVSPBCG_MAXL=5."""
	ret = cvodes.CVBPSpbcg(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpbcg() failed with flag %i"%(ret))
cvodes.CVBPSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBPSpbcg.restype = ctypes.c_int

def CVBPSpgmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpgmr links the CVBANDPPRE preconditioner to the CVSPGMR linear solver. It performs the following actions: 
	1) Calls the CVSPGMR specification routine and attaches the CVSPGMR linear solver to the integrator memory; 
	2) Sets the preconditioner data structure for CVSPGMR 
	3) Sets the preconditioner setup routine for CVSPGMR 
	4) Sets the preconditioner solve routine for CVSPGMR 

	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype	(int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in sundials_iterative.h.  These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPGMR solver. Pass 0 to use the default value CVSPGMR_MAXL=5."""
	ret = cvodes.CVBPSpgmr(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpgmr() failed with flag %i"%(ret))
cvodes.CVBPSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvodes.CVBPSpgmr.restype = ctypes.c_int

def CVBandPrecFree(bp_data):
 	"""CVBandPrecFree frees the memory allocated by CVBandPrecAlloc in the argument bp_data."""
	ret = cvodes.CVBandPrecFree(bp_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecFree() failed with flag %i"%(ret))
cvodes.CVBandPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
cvodes.CVBandPrecFree.restype = None

def CVBandPrecGetWorkSpace(bp_data):
	"""CVBandPrecGetWorkSpace returns the CVBandPrec real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvodes.CVBandPrecGetWorkSpace(bp_data, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVBandPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVBandPrecGetWorkSpace.restype = ctypes.c_int

def CVBandPrecGetNumRhsEvals(bp_data):
 	"""CVBandPrecGetNumGfnEvals returns the number of calls made from CVBANDPRE to the user's RHS."""
	nfevalsBP = ctypes.c_long(0)
	ret = cvodes.CVBandPrecGetNumRhsEvals(bp_data, ctypes.byref(nfevalsBP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecGetNumRhsEvals() failed with flag %i"%(ret))
	return nfevalsBP.value
cvodes.CVBandPrecGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVBandPrecGetNumRhsEvals.restype = ctypes.c_int

def CVBandPrecGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBANDPRE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVBandPrecGetReturnFlagName(flag)
cvodes.CVBandPrecGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVBandPrecGetReturnFlagName.restype = ctypes.c_char_p

def CVBandPrecAllocB(cvadj_mem, nB, muB, mlB):
	ret = cvodes.CVBandPrecAllocB(cvadj_mem, nB, muB, mlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecAllocB() failed with flag %i"%(ret))
cvodes.CVBandPrecAllocB.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvodes.CVBandPrecAllocB.restype = ctypes.c_int

def CVBPSptfqmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBPSptfqmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSptfqmrB() failed with flag %i"%(ret))
cvodes.CVBPSptfqmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBPSptfqmrB.restype = ctypes.c_int

def CVBPSpbcgB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBPSpbcgB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpbcgB() failed with flag %i"%(ret))
cvodes.CVBPSpbcgB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBPSpbcgB.restype = ctypes.c_int

def CVBPSpgmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVBPSpgmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpgmrB() failed with flag %i"%(ret))
cvodes.CVBPSpgmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVBPSpgmrB.restype = ctypes.c_int

def CVBandPrecFreeB(cvadj_mem):
	ret = cvodes.CVBandPrecFreeB(cvadj_mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecFreeB() failed with flag %i"%(ret))
cvodes.CVBandPrecFreeB.argtypes = [ctypes.c_void_p]
cvodes.CVBandPrecFreeB.restype = None

#################
# cvodes_band.h #
#################

CVB_MSBJ = 50
CVB_DGMAX = 0.2

CVBAND_SUCCESS = 0
CVBAND_MEM_NULL = -1
CVBAND_LMEM_NULL = -2
CVBAND_ILL_INPUT = -3
CVBAND_MEM_FAIL = -4

#Additional last_flag values

CVBAND_JACFUNC_UNRECVR = -5
CVBAND_JACFUNC_RECVR = -6

#Return values for adjoint module

CVBAND_ADJMEM_NULL = -101
CVBAND_LMEMB_NULL = -102

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
	init (int)	dimension of matrix to create; always square
	mu (int)	width of portion of band above diagonal
	ml (int)	width of portion of band below diagonal
	smu (int)	is the storage upper bandwidth, mu <= smu <= size-1.  The BandGBTRF routine writes the LU factors into the storage for A. The upper triangular factor U, however, may have an upper bandwidth as big as MIN(size-1,mu+ml) because of partial pivoting. The smu field holds the upper bandwidth allocated for A."""
		if type(init) == int:
			self.size = init
			self.mu = mu
			self.ml = ml
			if smu < mu:
				raise ValueError("smu must be greater than or equal to mu")
			self.smu = smu
			self.data = cvodes.BandAllocMat(init, mu, ml, smu)
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
		"""x.__repr__(i, y) <==> repr(x)"""
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
				cvodes.BandFreeMat(self.data)
			except:
				pass

def BandAllocMat(N, mu, ml, smu):
	"""Allocates memory for a banded matrix. Should not be called directly, rather instatiate a BandMat object.
	N (int)		dimension of matrix to create; always square
	mu (int)	width of portion of band above diagonal
	ml (int)	width of portion of band below diagonal
	smu (int)	is the storage upper bandwidth, mu <= smu <= size-1.  The BandGBTRF routine writes the LU factors into the storage for A. The upper triangular factor U, however, may have an upper bandwidth as big as MIN(size-1,mu+ml) because of partial pivoting. The smu field holds the upper bandwidth allocated for A."""
	return BandMat(N, mu, ml, smu)
cvodes.BandAllocMat.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvodes.BandAllocMat.restype = ctypes.POINTER(_BandMat)

def BandAllocPiv(N):
	"""BandAllocPiv allocates memory for pivot information to be filled in by the BandGBTRF routine during the factorization of an N by N band matrix. Returns a pointer, which should be passed as is to the [CV]BandPiv* functions.
	N (int)	the size of the banded matrix for which to allocate pivot information space [int]"""
	return cvodes.BandAllocPiv(N)
cvodes.BandAllocPiv.argtypes = [ctypes.c_int]
cvodes.BandAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def BandGBTRF(A, p):
	"""BandGBTRF performs the LU factorization of the N by N band matrix A. This is done using standard Gaussian elimination with partial pivoting.
	
	A (BandMat)	The banded matrix to factorize
	p (*long int)	The array of pivots

A successful LU factorization leaves the "matrix" A and the pivot array p with the following information:
	1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.
	2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower triangular matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower triangular part of A contains the multipliers, I-L.

BandGBTRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero.

Important Note: 'A' must be allocated to accommodate the increase in upper bandwidth that occurs during factorization. If mathematically, A is a band matrix with upper bandwidth mu and lower bandwidth ml, then the upper triangular factor U can have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower triangular factor L has lower bandwidth ml. Allocate A with call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are as defined above. The user does not have to zero the "extra" storage allocated for the purpose of factorization. This will handled by the BandGBTRF routine.  ret = cvodes.BandGBTRF(A, p)"""
	ret = cvodes.BandGBTRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRF() failed with flag %i"%(ret))
cvodes.BandGBTRF.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long)]
cvodes.BandGBTRF.restype = ctypes.c_long

def BandGBTRS(A, p, b):
	"""BandGBTRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in BandGBTRF. The solution x is returned in b. This routine cannot fail if the corresponding call to BandGBTRF did not fail.
	A (BandMat)	The banded matrix to factorize
	p (*long int)	The array of pivots
	b (*realtype)	The array containing the right hand side of the system A x = b"""

	ret = cvodes.BandGBTRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRS() failed with flag %i"%(ret))
cvodes.BandGBTRS.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
cvodes.BandGBTRS.restype = None

def BandZero(A):
	"""Zeroes the entire matrix
	A (BandMat)	the matrix"""
	cvodes.BandZero(A.data)
cvodes.BandZero.argtypes = [_BandMat] 
cvodes.BandZero.restype = None

def BandCopy(A, B, copymu, copyml):
	"""Copies the contents of banded matrix B into A, overwriting A's contents.
	A (BandMat)	the destination matrix
	B (BandMat)	the source matrix
	copymu (int)	upper band width of matrices
	copyml (int)	lower band width of matrices"""
	cvodes.BandCopy(A.data, B.data, copymu, copyml)
cvodes.BandCopy.argtypes = [_BandMat, _BandMat, ctypes.c_long, ctypes.c_long]
cvodes.BandCopy.restype = None

def BandScale(c, A):
	"""Scales the matrix A by c
	A (BandMat)	the matrix
	c (float)	the scale factor"""
	cvodes.BandScale(c, A.data)
cvodes.BandScale.argtypes = [realtype, _BandMat]
cvodes.BandScale.restype = None

def BandAddI(A):
	"""Adds 1.0 to the diagonal of the banded matrix A
	A (BandMat)	the matrix"""
	cvodes.BandAddI(A.data)
cvodes.BandAddI.argtypes = [_BandMat]
cvodes.BandAddI.restype = None

def BandFreeMat(A):
	"""BandFreeMat frees the memory allocated by BandAllocMat for the band matrix A.
	A (BandMat)	the matrix"""
	cvodes.BandFreeMat(A.data)
	del A.data
cvodes.BandFreeMat.argtypes = [_BandMat]
cvodes.BandFreeMat.restype = None

def BandFreePiv(p):
	"""BandFreePiv frees the pivot information storage memory p allocated by BandAllocPiv.
	p (*long int)	the pivot array"""
	cvodes.BandFreePiv(p)
cvodes.BandFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvodes.BandFreePiv.restype = None

def BandPrint(A):
	"""Print out the banded matix A
	A (BandMat)	the matrix"""
	cvodes.BandPrint(A.data)
cvodes.BandPrint
cvodes.BandPrint.restype = None

CVBandJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.POINTER(_BandMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVBandJacFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBandSetJac.

The callable python object must take exactly 11 parameters, which will be passed in as
	N (int)			the length of all NVector arguments
	mupper (int)		upper band width
	mlower (int)		lower band width
	J (BandMat)		the matrix that will be loaded with anpproximation of the Jacobian Matrix J = (df_i/dy_j) at the point (t,y).
	t (realtype)		the current value of the independent variable
	y (NVector)		the current value of the dependent variable vector
	fy (NVector)		f(t,y)
	jac_data (c_void_p)	pointer to user data set by CVBandSetJacFunc
	tmp1 (NVector)		preallocated temporary working space
	tmp2 (NVector)		preallocated temporary working space
	tmp3 (NVector)		preallocated temporary working space

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
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
	N (int)				the size of the ODE system.  
	mupper (int)			the upper bandwidth of the band Jacobian approximation.
	mlower (int)			the lower bandwidth of the band Jacobian approximation."""
	ret = cvodes.CVBand(cvodememobj.obj, N, mupper, mlower)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvodes.CVBand.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvodes.CVBand.restype = ctypes.c_int

def CVBandSetJacFn(cvodememobj, func, jac_data):
	"""CVBandSetJacFn sets the band Jacobian approximation function to be used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func	is a python callable taking the parameters
		N (int)			the length of all NVector arguments
		mupper (int)		upper band width
		mlower (int)		lower band width
		J (BandMat)		the matrix that will be loaded with anpproximation of the Jacobian Matrix J = (df_i/dy_j) at the point (t,y).
		t (realtype)		the current value of the independent variable
		y (NVector)		the current value of the dependent variable vector
		fy (NVector)		f(t,y)
		jac_data (c_void_p)	pointer to user data set by CVBandSetJacFunc
		tmp1 (NVector)		preallocated temporary working space
		tmp2 (NVector)		preallocated temporary working space
		tmp3 (NVector)		preallocated temporary working space
	jac_data			a pointer to user data"""
	ret = cvodes.CVBandSetJacFn(cvodememobj.obj, WrapCallbackCVBandJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvodes.CVBandSetJacFn.argtypes = [ctypes.c_void_p, CVBandJacFn, ctypes.c_void_p]
cvodes.CVBand.restype = ctypes.c_int

def CVBandGetWorkSpace(cvodememobj):
	"""CVBandGetWorkSpace returns the CVBand real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long()
	leniwLS = ctypes.c_long()
	ret = cvodes.CVBandGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVBandGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVBandGetWorkSpace.restype = ctypes.c_int

def CVBandGetNumJacEvals(cvodememobj):
 	"""CVBandGetNumGfnEvals returns the number of calls made from CVBAND to the user's jacobian function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVBandGetNumJacEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumJacEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVBandGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVBandGetNumJacEvals.restype = ctypes.c_int

def CVBandGetNumRhsEvals(cvodememobj):
 	"""CVBandGetNumGfnEvals returns the number of calls made from CVBAND to the user's RHS function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVBandGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVBandGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVBandGetNumRhsEvals.restype = ctypes.c_int

def CVBandGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVBAND interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int()
	ret = cvodes.CVBandGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVBandGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVBandGetLastFlag.restype = ctypes.c_int

def CVBandGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBAND return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVBandGetReturnFlagName(flag)
cvodes.CVBandGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVBandGetReturnFlagName.restype = ctypes.c_char_p

CVBandJacFnB = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.POINTER(_BandMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVBandJacFnB(func):
	exec 'def __CallbackInterface_%s(nB, mupperB, mlowerB, JB, t, y, yB, fyB, jac_dataB, tmp1B, tmp2B, tmp3B):\n\treturn __ActualCallback[%i](nB, mupperB, mlowerB, BandMat(JB), t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(fyB), jac_dataB, nvecserial.NVector(tmp1B), nvecserial.NVector(tmp2B), nvecserial.NVector(tmp3B))'%(func.func_name, len(__ActualCallback))
	if func == None:
		return ctypes.cast(None, CVBandJacFnB)
	__ActualCallback.append(func)
	tmp = CVBandJacFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVBandB(cvadj_mem, nB, mupperB, mlowerB):
	"""CVBandB links the main CVODES integrator with the CVBAND linear solver for the backward integration."""
	ret = cvodes.CVBandB(cvadj_mem, nB, mupperB, mlowerB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandB() failed with flag %i"%(ret))
cvodes.CVBandB.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvodes.CVBandB.restype = ctypes.c_int

def CVBandSetJacFnB(cvadj_mem, bjacB, jac_dataB):
	ret = cvodes.CVBandSetJacFnB(cvadj_mem, WrapCallbackCVBandJacFnB(bjacB), jac_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandSetJacFnB() failed with flag %i"%(ret))
cvodes.CVBandSetJacFnB.argtypes = [ctypes.c_void_p, CVBandJacFnB, ctypes.c_void_p]
cvodes.CVBandSetJacFnB.restype = ctypes.c_int

##################
# cvodes_dense.h #
##################

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
		"""x.__repr__(i, y) <==> repr(x)"""
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
	M (int)	number of rows
	N (int)	number of columns"""
		if type(M) == int and N is not None and type(N) == int:
			self.M = N
			self.N = M
			self.data = cvodes.DenseAllocMat(N,M)
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
		"""x.__repr__(i, y) <==> repr(x)"""
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
				cvodes.DenseFreeMat(self.data)
			except:
				pass

def DenseAllocMat(M, N):
	"""Allocates memory for a dense matrix. Should not be called directly, rather instatiate a DenseMat object.
	M (int)	number of rows
	N (int)	number of colums"""
	return DenseMat(M, N)
cvodes.DenseAllocMat.argtypes = [ctypes.c_int, ctypes.c_int]
cvodes.DenseAllocMat.restype = ctypes.POINTER(_DenseMat)

def DenseAllocPiv(N):
	"""DenseAllocPiv allocates memory for pivot information to be filled in by the DenseGETRF routine during the factorization of an N by N dense matrix. Returns a pointer, which should be passed as is to the [CV]DensePiv* functions.
	N (int)	the size of the dense matrix for which to allocate pivot information space [int]"""
	return cvodes.DenseAllocPiv(N)
cvodes.DenseAllocPiv.argtypes = [ctypes.c_int]
cvodes.DenseAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def DenseGETRF(A, p):
	"""DenseGETRF performs the LU factorization of the M by N dense matrix A. This is done using standard Gaussian elimination with partial (row) pivoting. Note that this only applies to matrices with M >= N and full column rank.

	A (DenseMat)	The matrix to factorize
	p (*long int)	The array of pivots

A successful LU factorization leaves the matrix A and the pivot array p with the following information:
	1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.
	2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower trapezoidal part of A contains the multipliers, I-L.

For square matrices (M=N), L is unit lower triangular.

DenseGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	ret = cvodes.DenseGETRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRF() failed with flag %i"%(ret))
cvodes.DenseGETRF.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long)]
cvodes.DenseGETRF.restype = ctypes.c_long

def DenseGETRS(A, p, b):
	"""DenseGETRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in DenseGETRF. The solution x is returned in b. This routine cannot fail if the corresponding call to DenseGETRF did not fail.
	A (DenseMat)	The banded matrix to factorize
	p (*long int)	The array of pivots
	b (*realtype)	The array containing the right hand side of the system A x = b

DenseGETRS does NOT check for a squre matrix!"""
	ret = cvodes.DenseGETRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRS() failed with flag %i"%(ret))
cvodes.DenseGETRS.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
cvodes.DenseGETRS.restype = None

def DenseZero(A):
	"""Zeroes the entire matrix
	A (DenseMat)	the matrix"""
	cvodes.DenseZero(A.data)
cvodes.DenseZero.argtypes = [_DenseMat] 
cvodes.DenseZero.restype = None

def DenseCopy(A, B):
	"""Copies the contents of dense matrix B into A, overwriting A's contents.
	A (DenseMat)	the destination matrix
	B (DenseMat)	the source matrix
	copymu (int)	upper band width of matrices
	copyml (int)	lower band width of matrices"""
	cvodes.DenseCopy(A.data, B.data)
cvodes.DenseCopy.argtypes = [_DenseMat, _DenseMat]
cvodes.DenseCopy.restype = None

def DenseScale(c, A):
	"""Scales the matrix A by c
	A (DenseMat)	the matrix
	c (float)	the scale factor"""
	cvodes.DenseScale(c, A.data)
cvodes.DenseScale.argtypes = [realtype, _DenseMat]
cvodes.DenseScale.restype = None

def DenseAddI(A):
	"""Adds 1.0 to the diagonal of the dense matrix A
	A (DenseMat)	the matrix"""
	cvodes.DenseAddI(A.data)
cvodes.DenseAddI.argtypes = [_DenseMat]
cvodes.DenseAddI.restype = None

def DenseFreeMat(A):
	"""DenseFreeMat frees the memory allocated by DenseAllocMat for the band matrix A.
	A (DenseMat)	the matrix"""
	cvodes.DenseFreeMat(A.data)
	del A.data
cvodes.DenseFreeMat.argtypes = [_DenseMat]
cvodes.DenseFreeMat.restype = None

def DenseFreePiv(p):
	"""DenseFreePiv frees the pivot information storage memory p allocated by DenseAllocPiv.
	p (*long int)	The array of pivots"""
	cvodes.DenseFreePiv(p)
cvodes.DenseFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvodes.DenseFreePiv.restype = None

def DensePrint(A):
	"""Print out the dense matix A
	A (DenseMat)	the matrix"""
	cvodes.DensePrint(A.data)
cvodes.DensePrint
cvodes.DensePrint.restype = None

CVDenseJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(_DenseMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVDenseJacFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBandSetJac.

The callable python object must take exactly 9 parameters, which will be passed in as
	N (int)			the length of all NVector arguments
	J (BandMat)		the matrix that will be loaded with anpproximation of the Jacobian Matrix J = (df_i/dy_j) at the point (t,y).
	t (realtype)		the current value of the independent variable
	y (NVector)		the current value of the dependent variable vector
	fy (NVector)		f(t,y)
	jac_data (c_void_p)	pointer to user data set by CVDenseSetJacFunc
	tmp1 (NVector)		preallocated temporary working space
	tmp2 (NVector)		preallocated temporary working space
	tmp3 (NVector)		preallocated temporary working space

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVDenseJacFn)
	exec 'def __CallbackInterface_%s(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](N, DenseMat(J), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVDenseJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVDense(cvodememobj, N):
	"""A call to the CVDense function links the main CVODE integrator with the CVDENSE linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	N (int)				the size of the ODE system."""
	ret = cvodes.CVDense(cvodememobj.obj, N)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvodes.CVDense.argtypes = [ctypes.c_void_p, ctypes.c_long]
cvodes.CVDense.restype = ctypes.c_int

def CVDenseSetJacFn(cvodememobj, func, jac_data):
	"""CVDenseSetJacFn sets the dense Jacobian approximation function to be used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	func	is a python callable taking the parameters
		N (int)			the length of all NVector arguments
		J (BandMat)		the matrix that will be loaded with anpproximation of the Jacobian Matrix J = (df_i/dy_j) at the point (t,y).
		t (realtype)		the current value of the independent variable
		y (NVector)		the current value of the dependent variable vector
		fy (NVector)		f(t,y)
		jac_data (c_void_p)	pointer to user data set by CVDenseSetJacFunc
		tmp1 (NVector)		preallocated temporary working space
		tmp2 (NVector)		preallocated temporary working space
		tmp3 (NVector)		preallocated temporary working space
	jac_data			a pointer to user data"""
	ret = cvodes.CVDenseSetJacFn(cvodememobj.obj, WrapCallbackCVDenseJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvodes.CVDenseSetJacFn.argtypes = [ctypes.c_void_p, CVDenseJacFn, ctypes.c_void_p]
cvodes.CVDense.restype = ctypes.c_int

def CVDenseGetWorkSpace(cvodememobj):
	"""CVDenseGetWorkSpace returns the CVDense real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long()
	leniwLS = ctypes.c_long()
	ret = cvodes.CVDenseGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVDenseGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVDenseGetWorkSpace.restype = ctypes.c_int

def CVDenseGetNumJacEvals(cvodememobj):
 	"""CVDenseGetNumGfnEvals returns the number of calls made from CVDENSE to the user's Jacobian function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVDenseGetNumJacEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumJacEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVDenseGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVDenseGetNumJacEvals.restype = ctypes.c_int

def CVDenseGetNumRhsEvals(cvodememobj):
 	"""CVDenseGetNumGfnEvals returns the number of calls made from CVDENSE to the user's RHS function.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvodes.CVDenseGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVDenseGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVDenseGetNumRhsEvals.restype = ctypes.c_int

def CVDenseGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVDENSE interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	retval = ctypes.c_int()
	ret = cvodes.CVDenseGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvodes.CVDenseGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVDenseGetLastFlag.restype = ctypes.c_int

def CVDenseGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDENSE return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVDenseGetReturnFlagName(flag)
cvodes.CVDenseGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVDenseGetReturnFlagName.restype = ctypes.c_char_p

CVDenseJacFnB = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(_DenseMat), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVDenseJacFnB(func):
	if func == None:
		return ctypes.cast(None, CVDenseJacFnB)
	exec 'def __CallbackInterface_%s(nB, JB, t, y, yB, fyB, jac_dataB, tmp1B, tmp2B, tmp3B):\n\treturn __ActualCallback[%i](nB, DenseMat(JB), t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(fyB), jac_dataB, nvecserial.NVector(tmp1B), nvecserial.NVector(tmp2B), nvecserial.NVector(tmp3B))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVDenseJacFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVDenseB(cvadj_mem, nB):
	"""CVDenseB links the main CVODES integrator with the CVDENSE linear solver for the backward integration."""
	ret = cvodes.CVDenseB(cvadj_mem, nB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseB() failed with flag %i"%(ret))
cvodes.CVDenseB.argtypes = [ctypes.c_void_p, ctypes.c_long]
cvodes.CVDenseB.restype = ctypes.c_int

def CVDenseSetJacFnB(cvadj_mem, djacB, jac_dataB):
	ret = cvodes.CVDenseSetJacFnB(cvadj_mem, WrapCallbackCVDenseJacFnB(djacB), jac_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseSetJacFnB() failed with flag %i"%(ret))
cvodes.CVDenseSetJacFnB.argtypes = [ctypes.c_void_p, CVDenseJacFnB, ctypes.c_void_p]
cvodes.CVDenseSetJacFnB.restype = ctypes.c_int

##################
# cvodes_spils.h #
##################

CVSPILS_MAXL = 5
CVSPILS_MSBPRE = 50
CVSPILS_DGMAX = 0.2
CVSPILS_DELT = 0.05

CVSPILS_SUCCESS = 0
CVSPILS_MEM_NULL = -1
CVSPILS_LMEM_NULL = -2
CVSPILS_ILL_INPUT = -3
CVSPILS_MEM_FAIL = -4

#Return values for the adjoint module
CVSPILS_ADJMEM_NULL = -101
CVSPILS_LMEMB_NULL = -102

CVSpilsPrecSetupFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, ctypes.POINTER(ctypes.c_int), realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSetupFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBandSetJac.

The callable python object must take exactly 10 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
	y (NVector)		the current value of the dependent variable vector
	fy (NVector)		f(t,y)
	jok (*int)		set to 1 or 0 to indicate whether the jacobian needs to be recomputed from scratch (0) or not (1)
	jCurPtr (*int)		should be set to 1 if Jacobian data was recomputed, otherwise 0
	gamma (realtype)	the scalar appearing in the Newton matrix
	P_data (c_void_p)	pointer to user data set by CVSpilsSetPreconditioner
	tmp1 (NVector)		preallocated temporary working space
	tmp2 (NVector)		preallocated temporary working space
	tmp3 (NVector)		preallocated temporary working space

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSetupFn)
	exec 'def __CallbackInterface_%s(t, y, fy, jok, jcurPtr, gamma, P_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(fy), jok, jcurPtr, gamma, P_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSetupFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsPrecSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, realtype, ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSolveFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBandSetJac.

The callable python object must take exactly 10 parameters, which will be passed in as
	t (realtype)		the current value of the independent variable
	y (NVector)		the current value of the dependent variable vector
	fy (NVector)		f(t,y)
	r (NVector)		right hand side vector of the linear system
	z (NVector)		the ouput vector computed by PrecSolve
	gamma (realtype)	the scalar appearing in the Newton matrix
	delta (realtype)	input tolerance for use by PSolve
	lr (int)		use left preconditioner (1) or right preconditioner (2)
	P_data (c_void_p)	pointer to user data set by CVSpilsSetPreconditioner
	tmp (NVector)		preallocated temporary working space

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSolveFn)
	exec 'def __CallbackInterface_%s(t, y, fy, r, z, gamma, delta, lr, P_data, tmp):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(fy), nvecserial.NVector(r), nvecserial.NVector(z), gamma, delta, lr, P_data, nvecserial.NVector(tmp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsJacTimesVecFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsJacTimesVecFn(func):
	"""Returns a callback wrapper around the given python callable object (func). This function should never be called directly, as it is called implicitly by CVBandSetJac.

The callable python object must take exactly 7 parameters, which will be passed in as
	v (NVector)		the NVector to be multiplied by J
	Jv			the output NVector containing J*v
	t (realtype)		the current value of the independent variable
	y (NVector)		the current value of the dependent variable vector
	fy (NVector)		f(t,y)
	jac_data (c_void_p)	pointer to user data set by CVSpilsSetJacTimesVecFn
	tmp (NVector)		preallocated temporary working space

and must return an integer of 0 in the case of no error, otherwise a user defined integer indicating an error condition."""
	if func == None:
		return ctypes.cast(None, CVSpilsJacTimesVecFn)
	exec 'def __CallbackInterface_%s(v, Jv, t, y, fy, jac_data, tmp):\n\treturn __ActualCallback[%i](nvecserial.NVector(v), nvecserial.NVector(Jv), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsJacTimesVecFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVSpilsSetPrecType(cvodememobj, pretype):
	"""CVSpilsSetPrecType (re)sets the type of preconditioner.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	pretype (int)			This must be one of PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH."""
	ret = cvodes.CVSpilsSetPrecType(cvodememobj.obj, pretype)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPrecType() failed with flag %i"%(ret))
cvodes.CVSpilsSetPrecType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetPrecType.restype = ctypes.c_int

def CVSpilsSetGSType(cvodememobj, gstype):
	"""CVSpilsSetGSType specifies the type of Gram-Schmidt orthogonalization to be used. 
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	gstype (int)			This must be one of the two enumeration constants MODIFIED_GS or CLASSICAL_GS defined in iterative.h. These correspond to using modified Gram-Schmidt and classical Gram-Schmidt, respectively.  Default value is MODIFIED_GS."""
	ret = cvodes.CVSpilsSetGSType(cvodememobj.obj, gstype)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetGSType() failed with flag %i"%(ret))
cvodes.CVSpilsSetGSType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetGSType.restype = ctypes.c_int

def CVSpilsSetMaxl(cvodememobj, maxl):
	"""CVSpilsSetMaxl (re)sets the maximum Krylov subspace size.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	maxl (int)			the maximum Krylov subspace size. A value <= 0, gives the default value."""
	ret = cvodes.CVSpilsSetMaxl(cvodememobj.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetMaxl() failed with flag %i"%(ret))
cvodes.CVSpilsSetMaxl.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetMaxl.restype = ctypes.c_int

def CVSpilsSetDelt(cvodememobj, delt):
	"""CVSpilsSetDelt specifies the factor by which the tolerance on the nonlinear iteration is multiplied to get a tolerance on the linear iteration.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()
	delt (realtype)			Default value is 0.05."""
	ret = cvodes.CVSpilsSetDelt(cvodememobj.obj, delt)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetDelt() failed with flag %i"%(ret))
cvodes.CVSpilsSetDelt.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVSpilsSetDelt.restype = ctypes.c_int

def CVSpilsSetPreconditioner(cvodememobj, pset, psolve, P_data):
	"""CVDenseSetJacFn sets the dense Jacobian approximation function to be used.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pset	is a python callable taking the parameters
		t (realtype)		the current value of the independent variable
		y (NVector)		the current value of the dependent variable vector
		fy (NVector)		f(t,y)
		jok (*int)		set to 1 or 0 to indicate whether the jacobian needs to be recomputed from scratch (0) or not (1)
		jCurPtr (*int)		should be set to 1 if Jacobian data was recomputed, otherwise 0
		gamma (realtype)	the scalar appearing in the Newton matrix
		P_data (c_void_p)	pointer to user data set by CVSpilsSetPreconditioner
		tmp1 (NVector)		preallocated temporary working space
		tmp2 (NVector)		preallocated temporary working space
		tmp3 (NVector)		preallocated temporary working space
	psolve	is a python callable taking the parameters
		t (realtype)		the current value of the independent variable
		y (NVector)		the current value of the dependent variable vector
		fy (NVector)		f(t,y)
		r (NVector)		right hand side vector of the linear system
		z (NVector)		the ouput vector computed by PrecSolve
		gamma (realtype)	the scalar appearing in the Newton matrix
		delta (realtype)	input tolerance for use by PSolve
		lr (int)		use left preconditioner (1) or right preconditioner (2)
		P_data (c_void_p)	pointer to user data set by CVSpilsSetPreconditioner
		tmp (NVector)		preallocated temporary working space
	P_data				a pointer to user data"""
	ret = cvodes.CVSpilsSetPreconditioner(cvodememobj.obj, WrapCallbackCVSpilsPrecSetupFn(pset), WrapCallbackCVSpilsPrecSolveFn(psolve), P_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPreconditioner() failed with flag %i"%(ret))
cvodes.CVSpilsSetPreconditioner.argtypes = [ctypes.c_void_p, CVSpilsPrecSetupFn, CVSpilsPrecSolveFn, ctypes.c_void_p]
cvodes.CVSpilsSetPreconditioner.restype = ctypes.c_int

def CVSpilsSetJacTimesVecFn(cvodememobj, jtimes, jac_data):
	"""CVSpilsSetJacTimesVecFn specifies the jtimes function and a pointer to user Jacobian data. This pointer is passed to jtimes every time the jtimes routine is called.  Default is to use an internal finite difference approximation routine.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	jtimes	is a python callable taking the parameters
		v (NVector)		the NVector to be multiplied by J
		Jv			the output NVector containing J*v
		t (realtype)		the current value of the independent variable
		y (NVector)		the current value of the dependent variable vector
		fy (NVector)		f(t,y)
		jac_data (c_void_p)	pointer to user data set by CVSpilsSetJacTimesVecFn
		tmp (NVector)		preallocated temporary working space
	jac_data			a pointer to user data"""
	ret = cvodes.CVSpilsSetJacTimesVecFn(cvodememobj.obj, WrapCallbackCVSpilsJacTimesVecFn(jtimes), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetJacTimesVecFn() failed with flag %i"%(ret))
cvodes.CVSpilsSetJacTimesVecFn.argtypes = [ctypes.c_void_p, CVSpilsJacTimesVecFn, ctypes.c_void_p]
cvodes.CVSpilsSetJacTimesVecFn.restype = ctypes.c_int

def CVSpilsGetWorkSpace(cvodememobj):
	"""CVSpilsGetWorkSpace returns the CVSpils real and integer workspaces as a tuple (integer, real)
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvodes.CVSpilsGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetWorkSpace.restype = ctypes.c_int

def CVSpilsGetNumPrecEvals(cvodememobj):
	"""CVSpilsGetNumPrecEvals returns the number of preconditioner evaluations, i.e. the number of calls made to PrecSetup with jok==FALSE.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	npevals = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumPrecEvals(cvodememobj.obj, ctypes.byref(npevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecEvals() failed with flag %i"%(ret))
	return npevals.value
cvodes.CVSpilsGetNumPrecEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumPrecEvals.restype = ctypes.c_int

def CVSpilsGetNumPrecSolves(cvodememobj):
	"""CVSpilsGetNumPrecSolves returns the number of calls made to PrecSolve.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	npsolves = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumPrecSolves(cvodememobj.obj, ctypes.byref(npsolves))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecSolves() failed with flag %i"%(ret))
	return npsolves.value
cvodes.CVSpilsGetNumPrecSolves.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumPrecSolves.restype = ctypes.c_int

def CVSpilsGetNumLinIters(cvodememobj):
	"""CVSpilsGetNumLinIters returns the number of linear iterations.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nliters = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumLinIters(cvodememobj.obj, ctypes.byref(nliters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumLinIters() failed with flag %i"%(ret))
	return nliters.value
cvodes.CVSpilsGetNumLinIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumLinIters.restype = ctypes.c_int

def CVSpilsGetNumConvFails(cvodememobj):
	"""CVSpilsGetNumConvFails returns the number of linear convergence failures.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nlcfails = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumConvFails(cvodememobj.obj, ctypes.byref(nlcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumConvFails() failed with flag %i"%(ret))
	return nlcfails.value
cvodes.CVSpilsGetNumConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumConvFails.restype = ctypes.c_int

def CVSpilsGetNumJtimesEvals(cvodememobj):
	"""CVSpilsGetNumJtimesEvals returns the number of calls to jtimes.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	njvevals = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumJtimesEvals(cvodememobj.obj, ctypes.byref(njvevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumJtimesEvals() failed with flag %i"%(ret))
	return njvevals.value
cvodes.CVSpilsGetNumJtimesEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumJtimesEvals.restype = ctypes.c_int

def CVSpilsGetNumRhsEvals(cvodememobj):
	"""CVSpilsGetNumRhsEvals returns the number of calls to the user f routine due to finite difference Jacobian times vector evaluation.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()"""
	nfevalsLS = ctypes.c_long(0)
	ret = cvodes.CVSpilsGetNumRhsEvals(cvodememobj.obj, ctypes.byref(nfevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumRhsEvals() failed with flag %i"%(ret))
	return nfevalsLS.value
cvodes.CVSpilsGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvodes.CVSpilsGetNumRhsEvals.restype = ctypes.c_int

def CVSpilsGetLastFlag(cvodememobj):
	"""Returns the last error flag set by any of the CVSPILS interface functions.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CVodeCreate()"""
	flag = ctypes.c_int(0)
	ret = cvodes.CVSpilsGetLastFlag(cvodememobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetLastFlag() failed with flag %i"%(ret))
	return flag.value
cvodes.CVSpilsGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvodes.CVSpilsGetLastFlag.restype = ctypes.c_int

def CVSpilsGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVSPILS return flag.
	flag (int)	the integer constant of which to get the name"""
	return cvodes.CVSpilsGetReturnFlagName(flag)
cvodes.CVSpilsGetReturnFlagName.argtypes = [ctypes.c_int]
cvodes.CVSpilsGetReturnFlagName.restype = ctypes.c_char_p

CVSpilsPrecSetupFnB = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, ctypes.POINTER(ctypes.c_int), realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSetupFnB(func):
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSetupFnB)
	exec 'def __CallbackInterface_%s(t, y, yB, fyB, jokB, jcurPtrB, gammaB, P_dataB, tmp1B, tmp2B, tmp3B):\nreturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(fyB), jokB, ctypes.byref(jcurPtrB), gammaB, P_dataB, nvecserial.NVector(tmp1B), nvecserial.NVector(tmp2B), nvecserial.NVector(tmp3B))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSetupFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsPrecSolveFnB = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, realtype, ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsPrecSolveFnB(func):
	if func == None:
		return ctypes.cast(None, CVSpilsPrecSolveFnB)
	exec 'def __CallbackInterface_%s(t, y, yB, fyB, rB, zB, gammaB, deltaB, lrB, P_dataB, tmpB):\nreturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(fyB), nvecserial.NVector(rB), nvecserial.NVector(zB), gammaB, deltaB, lrB, P_dataB, nvecserial.NVector(tmpB))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsPrecSolveFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVSpilsJacTimesVecFnB = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackCVSpilsJacTimesVecFnB(func):
	if func == None:
		return ctypes.cast(None, CVSpilsJacTimesVecFnB)
	exec 'def __CallbackInterface_%s(vB, JvB, t, y, yB, fyB, jac_dataB, tmpB):\nreturn __ActualCallback[%i](nvecserial.NVector(vB), nvecserial.NVector(JvB), t, nvecserial.NVector(y), nvecserial.NVector(yB), nvecserial.NVector(fyB), jac_dataB, nvecserial.NVector(tmpB))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVSpilsJacTimesVecFnB(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVSpilsSetPrecTypeB(cvadj_mem, pretypeB):
	ret = cvodes.CVSpilsSetPrecTypeB(cvadj_mem, pretypeB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPrecTypeB() failed with flag %i"%(ret))
cvodes.CVSpilsSetPrecTypeB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetPrecTypeB.restype = ctypes.c_int

def CVSpilsSetGSTypeB(cvadj_mem, gstypeB):
	ret = cvodes.CVSpilsSetGSTypeB(cvadj_mem, gstypeB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetGSTypeB() failed with flag %i"%(ret))
cvodes.CVSpilsSetGSTypeB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetGSTypeB.restype = ctypes.c_int

def CVSpilsSetDeltB(cvadj_mem, deltB):
	ret = cvodes.CVSpilsSetDeltB(cvadj_mem, deltB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetDeltB() failed with flag %i"%(ret))
cvodes.CVSpilsSetDeltB.argtypes = [ctypes.c_void_p, realtype]
cvodes.CVSpilsSetDeltB.restype = ctypes.c_int

def CVSpilsSetMaxlB(cvadj_mem, maxlB):
	ret = cvodes.CVSpilsSetMaxlB(cvadj_mem, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetMaxlB() failed with flag %i"%(ret))
cvodes.CVSpilsSetMaxlB.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvodes.CVSpilsSetMaxlB.restype = ctypes.c_int

def CVSpilsSetPreconditionerB(cvadj_mem, psetB, psolveB, P_dataB):
	ret = cvodes.CVSpilsSetPreconditionerB(cvadj_mem, WrapCallbackCVSpilsPrecSetupFnB(psetB), WrapCallbackCVSpilsPrecSolveFnB(psolveB), P_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetPreconditionerB() failed with flag %i"%(ret))
cvodes.CVSpilsSetPreconditionerB.argtypes = [ctypes.c_void_p, CVSpilsPrecSetupFnB, CVSpilsPrecSolveFnB, ctypes.c_void_p]
cvodes.CVSpilsSetPreconditionerB.restype = ctypes.c_int

def CVSpilsSetJacTimesVecFnB(cvadj_mem, jtimesB, jac_dataB):
	ret = cvodes.CVSpilsSetJacTimesVecFnB(cvadj_mem, WrapCallbackCVSpilsJacTimesVecFnB(jtimesB), jac_dataB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsSetJacTimesVecFnB() failed with flag %i"%(ret))
cvodes.CVSpilsSetJacTimesVecFnB.argtypes = [ctypes.c_void_p, CVSpilsJacTimesVecFnB, ctypes.c_void_p]
cvodes.CVSpilsSetJacTimesVecFnB.restype = ctypes.c_int

###################
# cvodes_spbcgs.h #
###################

def CVSpbcg(cvodememobj, pretype, maxl):
	"""A call to the CVSpbcg function links the main CVODE integrator with the CVSPBCG linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in iterative.h. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPBCG solver. Pass 0 to use the default value CVSPBCG_MAXL=5."""
	ret = cvodes.CVSpbcg(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpbcg() failed with flag %i"%(ret))
cvodes.CVSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSpbcg.restype = ctypes.c_int

def CVSpbcgB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVSpbcgB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpbcgB() failed with flag %i"%(ret))
cvodes.CVSpbcgB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSpbcgB.restype = ctypes.c_int

##################
# cvodes_spgmr.h #
##################

def CVSpgmr(cvodememobj, pretype, maxl):
	"""A call to the CVSpgmr function links the main CVODE integrator with the CVSPGMR linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype	(int)			the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in sundials_iterative.h.  These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPGMR solver. Pass 0 to use the default value CVSPGMR_MAXL=5."""
	ret = cvodes.CVSpgmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpgmr() failed with flag %i"%(ret))
cvodes.CVSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSpgmr.restype = ctypes.c_int

def CVSpgmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVSpgmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpgmrB() failed with flag %i"%(ret))
cvodes.CVSpgmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSpgmrB.restype = ctypes.c_int

####################
# cvodes_sptfqmr.h #
####################

def CVSptfqmr(cvodememobj, pretype, maxl):
	"""A call to the CVSptfqmr function links the main CVODE integrator with the CVSPTFQMR linear solver.
	cvodememobj (CVodeMemObj)	a CVodeMemObj as returned by CvodeCreate()
	pretype (int)			the type of user preconditioning to be done. This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.
	maxl (int)			the maximum Krylov dimension. This is an optional input to the CVSPTFQMR solver. Pass 0 to use the default value CVSPILS_MAXL=5."""
	ret = cvodes.CVSptfqmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSptfqmr() failed with flag %i"%(ret))
cvodes.CVSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSptfqmr.restype = ctypes.c_int

def CVSptfqmrB(cvadj_mem, pretypeB, maxlB):
	ret = cvodes.CVSptfqmrB(cvadj_mem, pretypeB, maxlB)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSptfqmrB() failed with flag %i"%(ret))
cvodes.CVSptfqmrB.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvodes.CVSptfqmrB.restype = ctypes.c_int
