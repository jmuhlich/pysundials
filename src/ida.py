#ida.py is part of the PySUNDIALS package, and is released under the
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
"""Python bindings for the ida, ida_band, ida_bbdpre, ida_dense, ida_spbcgs, ida_spmgr, ida_spils, and ida_sptfqmr header files."""
import ctypes
import sundials_core
import nvecserial

realtype = nvecserial.realtype
NVector = nvecserial.NVector

ida = sundials_core.loadlib("ida")

#Inputs to IDAMalloc, IDAReInit, IDACalcIC, and IDASolve.

#Relative and Absolute Tolerance Types
IDA_SS = 1
IDA_SV = 2
IDA_WF = 3

#Internal task job type
IDA_NORMAL = 1
IDA_ONE_STEP = 2
IDA_NORMAL_TSTOP = 3
IDA_ONE_STEP_TSTOP = 4

#icopt
IDA_YA_YDP_INIT = 1
IDA_Y_INIT = 2

#IDA return flags 
IDA_SUCCESS = 0
IDA_TSTOP_RETURN = 1
IDA_ROOT_RETURN = 2

IDA_WARNING = 99

IDA_MEM_NULL = -1
IDA_ILL_INPUT = -2
IDA_NO_MALLOC = -3
IDA_TOO_MUCH_WORK = -4
IDA_TOO_MUCH_ACC = -5
IDA_ERR_FAIL = -6
IDA_CONV_FAIL = -7
IDA_LINIT_FAIL = -8
IDA_LSETUP_FAIL = -9
IDA_LSOLVE_FAIL = -10
IDA_RES_FAIL = -11
IDA_CONSTR_FAIL = -12
IDA_REP_RES_ERR = -13

IDA_MEM_FAIL = -14

IDA_BAD_T = -15

IDA_BAD_EWT = -16
IDA_FIRST_RES_FAIL = -17
IDA_LINESEARCH_FAIL = -18
IDA_NO_RECOVERY = -19

IDA_RTFUNC_FAIL = -20

__Callback = []
__ActualCallback = []

class IdaMemObj(object):
	def __init__(self, obj):
		self.obj = obj

	def __del__(self):
		p = ctypes.c_void_p()
		p.value = self.obj
		ida.IDAFree(ctypes.byref(p))
ida.IDAFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
ida.IDAFree.restype = None

IDAResFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackIDAResFn(func):
	if func == None:
		return ctypes.cast(None, IDAResFn)
	exec 'def __CallbackInterface_%s(tt, yy, yp, rr, res_data):\n\treturn __ActualCallback[%i](tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), res_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDAResFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDARootFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype ), ctypes.c_void_p)
def WrapCallbackIDARootFn(func):
	if func == None:
		return ctypes.cast(None, IDARootFn)
	exec 'def __CallbackInterface_%s(t, y, yp, gout, g_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(yp), gout, g_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDARootFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDAEwtFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackIDAEwtFn(func):
	if func == None:
		return ctypes.cast(None, IDAEwtFn)
	exec 'def __CallbackInterface_%s(y, ewt, e_data):\n\treturn __ActualCallback[%i](nvecserial.NVector(y), nvecserial.NVector(ewt), e_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDAEwtFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDAErrHandlerFn = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_void_p)
def WrapCallbackIDAErrHandlerFn(func):
	if func == None:
		return ctypes.cast(None, IDAErrHandlerFn)
	exec 'def __CallbackInterface_%s(error_code, module, function, msg, eh_data):\n\treturn __ActualCallback[%i](error_code, module, function, msg, eh_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDAErrHandlerFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def IDACreate():
	obj = ida.IDACreate()
	if obj == None:
		raise AssertionError("SUNDIALS ERROR: IDACreate() failed - returned NULL pointer")
	return IdaMemObj(obj)
ida.IDACreate.argtypes = []
ida.IDACreate.restype = ctypes.c_void_p

def IDASetErrHandlerFn(ida_mem, ehfun, eh_data):
	ret = ida.IDASetErrHandlerFn(ida_mem.obj, WrapCallbackIDAErrHandlerFn(ehfun), eh_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetErrHandlerFn() failed with flag %i"%(ret))
ida.IDASetErrHandlerFn.argtypes = [ctypes.c_void_p, IDAErrHandlerFn, ctypes.c_void_p]
ida.IDASetErrHandlerFn.restype = ctypes.c_int

def IDASetErrFile(ida_mem, errfp):
	ret = ida.IDASetErrFile(ida_mem.obj, sundials_core.fdopen(errfp.fileno, errfp.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetErrFile() failed with flag %i"%(ret))
ida.IDASetErrFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
ida.IDASetErrFile.restype = ctypes.c_int

def IDASetRdata(ida_mem, res_data):
	ret = ida.IDASetRdata(ida_mem.obj, res_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetRdata() failed with flag %i"%(ret))
ida.IDASetRdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
ida.IDASetRdata.restype = ctypes.c_int

def IDASetEwtFn(ida_mem, efun, edata):
	ret = ida.IDASetEwtFn(ida_mem.obj, WrapCallbackIDAEwtFn(efun), edata)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetEwtFn() failed with flag %i"%(ret))
ida.IDASetEwtFn.argtypes = [ctypes.c_void_p, IDAEwtFn, ctypes.c_void_p]
ida.IDASetEwtFn.restype = ctypes.c_int

def IDASetMaxOrd(ida_mem, maxord):
	ret = ida.IDASetMaxOrd(ida_mem.obj, maxord)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxOrd() failed with flag %i"%(ret))
ida.IDASetMaxOrd.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxOrd.restype = ctypes.c_int

def IDASetMaxNumSteps(ida_mem, mxsteps):
	ret = ida.IDASetMaxNumSteps(ida_mem.obj, mxsteps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxNumSteps() failed with flag %i"%(ret))
ida.IDASetMaxNumSteps.argtypes = [ctypes.c_void_p, ctypes.c_long]
ida.IDASetMaxNumSteps.restype = ctypes.c_int

def IDASetInitStep(ida_mem, hin):
	ret = ida.IDASetInitStep(ida_mem.obj, hin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetInitStep() failed with flag %i"%(ret))
ida.IDASetInitStep.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetInitStep.restype = ctypes.c_int

def IDASetMaxStep(ida_mem, hmax):
	ret = ida.IDASetMaxStep(ida_mem.obj, hmax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxStep() failed with flag %i"%(ret))
ida.IDASetMaxStep.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetMaxStep.restype = ctypes.c_int

def IDASetStopTime(ida_mem, tstop):
	ret = ida.IDASetStopTime(ida_mem.obj, tstop)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetStopTime() failed with flag %i"%(ret))
ida.IDASetStopTime.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetStopTime.restype = ctypes.c_int

def IDASetNonlinConvCoef(ida_mem, epcon):
	ret = ida.IDASetNonlinConvCoef(ida_mem.obj, epcon)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetNonlinConvCoef() failed with flag %i"%(ret))
ida.IDASetNonlinConvCoef.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetNonlinConvCoef.restype = ctypes.c_int

def IDASetMaxErrTestFails(ida_mem, maxnef):
	ret = ida.IDASetMaxErrTestFails(ida_mem.obj, maxnef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxErrTestFails() failed with flag %i"%(ret))
ida.IDASetMaxErrTestFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxErrTestFails.restype = ctypes.c_int

def IDASetMaxNonlinIters(ida_mem, maxcor):
	ret = ida.IDASetMaxNonlinIters(ida_mem.obj, maxcor)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxNonlinIters() failed with flag %i"%(ret))
ida.IDASetMaxNonlinIters.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxNonlinIters.restype = ctypes.c_int

def IDASetMaxConvFails(ida_mem, maxncf):
	ret = ida.IDASetMaxConvFails(ida_mem.obj, maxncf)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxConvFails() failed with flag %i"%(ret))
ida.IDASetMaxConvFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxConvFails.restype = ctypes.c_int

def IDASetSuppressAlg(ida_mem, suppressalg):
	ret = ida.IDASetSuppressAlg(ida_mem.obj, ctypes.byref(suppressalg))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetSuppressAlg() failed with flag %i"%(ret))
ida.IDASetSuppressAlg.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetSuppressAlg.restype = ctypes.c_int

def IDASetId(ida_mem, id):
	ret = ida.IDASetId(ida_mem.obj, id.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetId() failed with flag %i"%(ret))
ida.IDASetId.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
ida.IDASetId.restype = ctypes.c_int

def IDASetConstraints(ida_mem, constraints):
	ret = ida.IDASetConstraints(ida_mem.obj, constraints.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetConstraints() failed with flag %i"%(ret))
ida.IDASetConstraints.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
ida.IDASetConstraints.restype = ctypes.c_int

def IDASetTolerances(ida_mem, itol, rtol, atol):
	if itol == IDA_SS:
		if type(atol) == realtype:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, ctypes.byref(atol))
		elif type(atol) == float or type(atol) == int:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, ctype.byref(realtype(atol)))
		elif atol == None:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, atol)
		else:
			raise TypeError("atol must be a floating point number if itol is IDA_SS")
	elif itol == IDA_SV:
		if type(atol) == NVector:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, atol.data)
		elif atol == None:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, atol)
		else:
			raise TypeError("atol must be an NVector if itol is IDA_SV")
	elif itol == IDA_WF:
		if atol == None:
			ret = ida.IDASetTolerances(ida_mem.obj, itol, rtol, atol)
		else:
			raise TypeError("atol must be None if itol is IDA_WF")
	else:
		raise ValueError("itol must be one of IDA_SS or IDA_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetTolerances() failed with flag %i"%(ret))
ida.IDASetTolerances.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype, ctypes.c_void_p]
ida.IDASetTolerances.restype = ctypes.c_int

def IDAMalloc(ida_mem, res, t0, yy0, yp0, itol, rtol, atol):
	if itol == IDA_SS:
		if type(atol) == realtype:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, ctypes.byref(atol))
		elif type(atol) == float or type(atol) == int:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, ctypes.byref(realtype(atol)))
		elif atol == None:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be a floating point number if itol is IDA_SS")
	elif itol == IDA_SV:
		if type(atol) == NVector:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol.data)
		elif atol == None:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be an NVector if itol is IDA_SV")
	elif itol == IDA_WF:
		if atol == None:
			ret = ida.IDAMalloc(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be None if itol is IDA_WF")
	else:
		raise ValueError("itol must be one of IDA_SS or IDA_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAMalloc() failed with flag %i"%(ret))
ida.IDAMalloc.argtypes = [ctypes.c_void_p, IDAResFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
ida.IDAMalloc.restype = ctypes.c_int

def IDAReInit(ida_mem, res, t0, yy0, yp0, itol, rtol, atol):
	if itol == IDA_SS:
		if type(atol) == realtype:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, ctypes.byref(atol))
		elif type(atol) == float or type(atol) == int:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, ctypes.byref(realtype(atol)))
		elif atol == None:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be a floating point number if itol is IDA_SS")
	elif itol == IDA_SV:
		if type(atol) == NVector:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol.data)
		elif atol == None:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be an NVector if itol is IDA_SV")
	elif itol == IDA_WF:
		if atol == None:
			ret = ida.IDAReInit(ida_mem.obj, WrapCallbackIDAResFn(res), t0, yy0.data, yp0.data, itol, rtol, atol)
		else:
			raise TypeError("atol must be None if itol is IDA_WF")
	else:
		raise ValueError("itol must be one of IDA_SS or IDA_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAReInit() failed with flag %i"%(ret))
ida.IDAReInit.argtypes = [ctypes.c_void_p, IDAResFn, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p]
ida.IDAReInit.restype = ctypes.c_int

def IDASetNonlinConvCoefIC(ida_mem, epiccon):
	ret = ida.IDASetNonlinConvCoefIC(ida_mem.obj, epiccon)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetNonlinConvCoefIC() failed with flag %i"%(ret))
ida.IDASetNonlinConvCoefIC.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetNonlinConvCoefIC.restype = ctypes.c_int

def IDASetMaxNumStepsIC(ida_mem, maxnh):
	ret = ida.IDASetMaxNumStepsIC(ida_mem.obj, maxnh)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxNumStepsIC() failed with flag %i"%(ret))
ida.IDASetMaxNumStepsIC.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxNumStepsIC.restype = ctypes.c_int

def IDASetMaxNumJacsIC(ida_mem, maxnj):
	ret = ida.IDASetMaxNumJacsIC(ida_mem.obj, maxnj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxNumJacsIC() failed with flag %i"%(ret))
ida.IDASetMaxNumJacsIC.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxNumJacsIC.restype = ctypes.c_int

def IDASetMaxNumItersIC(ida_mem, maxnit):
	ret = ida.IDASetMaxNumItersIC(ida_mem.obj, maxnit)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetMaxNumItersIC() failed with flag %i"%(ret))
ida.IDASetMaxNumItersIC.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetMaxNumItersIC.restype = ctypes.c_int

def IDASetLineSearchOffIC(ida_mem, lsoff):
	ret = ida.IDASetLineSearchOffIC(ida_mem.obj, ctypes.byref(lsoff))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetLineSearchOffIC() failed with flag %i"%(ret))
ida.IDASetLineSearchOffIC.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASetLineSearchOffIC.restype = ctypes.c_int

def IDASetStepToleranceIC(ida_mem, steptol):
	ret = ida.IDASetStepToleranceIC(ida_mem.obj, steptol)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASetStepToleranceIC() failed with flag %i"%(ret))
ida.IDASetStepToleranceIC.argtypes = [ctypes.c_void_p, realtype]
ida.IDASetStepToleranceIC.restype = ctypes.c_int

def IDARootInit(ida_mem, nrtfn, g, g_data):
	ret = ida.IDARootInit(ida_mem.obj, nrtfn, WrapCallbackIDARootFn(g), g_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDARootInit() failed with flag %i"%(ret))
ida.IDARootInit.argtypes = [ctypes.c_void_p, ctypes.c_int, IDARootFn, ctypes.c_void_p]
ida.IDARootInit.restype = ctypes.c_int

def IDACalcIC(ida_mem, icopt, tout1):
	ret = ida.IDACalcIC(ida_mem.obj, icopt, tout1)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDACalcIC() failed with flag %i"%(ret))
ida.IDACalcIC.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype]
ida.IDACalcIC.restype = ctypes.c_int

def IDASolve(ida_mem, tout, tret, yret, ypret, itask):
	ret = ida.IDASolve(ida_mem.obj, tout, tret, yret.data, ypret.data, itask)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASolve() failed with flag %i"%(ret))
	return ret
ida.IDASolve.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(realtype ), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int]
ida.IDASolve.restype = ctypes.c_int

def IDAGetSolution(ida_mem, t, yret, ypret):
	ret = ida.IDAGetSolution(ida_mem.obj, t, yret.data, ypret.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetSolution() failed with flag %i"%(ret))
ida.IDAGetSolution.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector)]
ida.IDAGetSolution.restype = ctypes.c_int

def IDAGetWorkSpace(ida_mem):
	lenrw = ctypes.c_long(0)
	leniw = ctypes.c_long(0)
	ret = ida.IDAGetWorkSpace(ida_mem.obj, ctypes.byref(lenrw), ctypes.byref(leniw))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetWorkSpace() failed with flag %i"%(ret))
	return (lenrw.value, leniw.value)
ida.IDAGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDAGetWorkSpace.restype = ctypes.c_int

def IDAGetNumSteps(ida_mem):
	nsteps = ctypes.c_long(0)
	ret = ida.IDAGetNumSteps(ida_mem.obj, ctypes.byref(nsteps))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumSteps() failed with flag %i"%(ret))
	return nsteps.value
ida.IDAGetNumSteps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumSteps.restype = ctypes.c_int

def IDAGetNumResEvals(ida_mem):
	nrevals = ctypes.c_long(0)
	ret = ida.IDAGetNumResEvals(ida_mem.obj, ctypes.byref(nrevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumResEvals() failed with flag %i"%(ret))
	return nrevals.value
ida.IDAGetNumResEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumResEvals.restype = ctypes.c_int

def IDAGetNumLinSolvSetups(ida_mem):
	nlinsetups = ctypes.c_long(0)
	ret = ida.IDAGetNumLinSolvSetups(ida_mem.obj, ctypes.byref(nlinsetups))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumLinSolvSetups() failed with flag %i"%(ret))
	return nlinsetups.value
ida.IDAGetNumLinSolvSetups.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumLinSolvSetups.restype = ctypes.c_int

def IDAGetNumErrTestFails(ida_mem):
	netfails = ctypes.c_long(0)
	ret = ida.IDAGetNumErrTestFails(ida_mem.obj, ctypes.byref(netfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumErrTestFails() failed with flag %i"%(ret))
	return netfails.value
ida.IDAGetNumErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumErrTestFails.restype = ctypes.c_int

def IDAGetNumBacktrackOps(ida_mem):
	nbacktr = ctypes.c_long(0)
	ret = ida.IDAGetNumBacktrackOps(ida_mem.obj, ctypes.byref(nbacktr))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumBacktrackOps() failed with flag %i"%(ret))
	return nbacktr.value
ida.IDAGetNumBacktrackOps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumBacktrackOps.restype = ctypes.c_int

def IDAGetConsistentIC(ida_mem, yy0_mod, yp0_mod):
	ret = ida.IDAGetConsistentIC(ida_mem.obj, yy0_mod.data, yp0_mod.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetConsistentIC() failed with flag %i"%(ret))
ida.IDAGetConsistentIC.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector)]
ida.IDAGetConsistentIC.restype = ctypes.c_int

def IDAGetLastOrder(ida_mem):
	klast = ctypes.c_int(0)
	ret = ida.IDAGetLastOrder(ida_mem.obj, ctypes.byref(klast))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetLastOrder() failed with flag %i"%(ret))
	return klast.value
ida.IDAGetLastOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDAGetLastOrder.restype = ctypes.c_int

def IDAGetCurrentOrder(ida_mem):
	kcur = ctypes.c_int(0)
	ret = ida.IDAGetCurrentOrder(ida_mem.obj, ctypes.byref(kcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetCurrentOrder() failed with flag %i"%(ret))
	return kcur.value
ida.IDAGetCurrentOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDAGetCurrentOrder.restype = ctypes.c_int

def IDAGetActualInitStep(ida_mem):
	hinused = ctypes.c_int(0)
	ret = ida.IDAGetActualInitStep(ida_mem.obj, ctypes.byref(hinused))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetActualInitStep() failed with flag %i"%(ret))
	return hinused.value
ida.IDAGetActualInitStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype )]
ida.IDAGetActualInitStep.restype = ctypes.c_int

def IDAGetLastStep(ida_mem):
	hlast = ctypes.c_double(0.0)
	ret = ida.IDAGetLastStep(ida_mem.obj, ctypes.byref(hlast))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetLastStep() failed with flag %i"%(ret))
	return hlast.value
ida.IDAGetLastStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype )]
ida.IDAGetLastStep.restype = ctypes.c_int

def IDAGetCurrentStep(ida_mem):
	hcur = ctypes.c_double(0.0)
	ret = ida.IDAGetCurrentStep(ida_mem.obj, ctypes.byref(hcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetCurrentStep() failed with flag %i"%(ret))
	return hcur.value
ida.IDAGetCurrentStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype )]
ida.IDAGetCurrentStep.restype = ctypes.c_int

def IDAGetCurrentTime(ida_mem):
	tcur = ctypes.c_double(0.0)
	ret = ida.IDAGetCurrentTime(ida_mem.obj, ctypes.byref(tcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetCurrentTime() failed with flag %i"%(ret))
	return tcur.value
ida.IDAGetCurrentTime.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype )]
ida.IDAGetCurrentTime.restype = ctypes.c_int

def IDAGetTolScaleFactor(ida_mem):
	tolsfact = ctypes.c_double(0.0)
	ret = ida.IDAGetTolScaleFactor(ida_mem.obj, ctypes.byref(tolsfact))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetTolScaleFactor() failed with flag %i"%(ret))
	return tolsfact.value
ida.IDAGetTolScaleFactor.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype )]
ida.IDAGetTolScaleFactor.restype = ctypes.c_int

def IDAGetErrWeights(ida_mem, eweight):
	ret = ida.IDAGetErrWeights(ida_mem.obj, eweight.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetErrWeights() failed with flag %i"%(ret))
ida.IDAGetErrWeights.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
ida.IDAGetErrWeights.restype = ctypes.c_int

def IDAGetEstLocalErrors(ida_mem, ele):
	ret = ida.IDAGetEstLocalErrors(ida_mem.obj, ele.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetEstLocalErrors() failed with flag %i"%(ret))
ida.IDAGetEstLocalErrors.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
ida.IDAGetEstLocalErrors.restype = ctypes.c_int

def IDAGetNumGEvals(ida_mem):
	ngevals = ctypes.c_long(0)
	ret = ida.IDAGetNumGEvals(ida_mem.obj, ctypes.byref(ngevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumGEvals() failed with flag %i"%(ret))
	return ngevals.value
ida.IDAGetNumGEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumGEvals.restype = ctypes.c_int

def IDAGetRootInfo(ida_mem, numroots):
	rootsfound = (ctypes.c_int*numroots)()
	ret = ida.IDAGetRootInfo(ida_mem.obj, rootsfound)
	return [rootsfound[i] for i in range(numroots)]
ida.IDAGetRootInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDAGetRootInfo.restype = ctypes.c_int

def IDAGetIntegratorStats(ida_mem):
	nsteps = ctypes.c_long(0)
	nrevals = ctypes.c_long(0)
	nlinsetups = ctypes.c_long(0)
	netfails = ctypes.c_long(0)
	qlast = ctypes.c_int(0)
	qcur = ctypes.c_int(0)
	hinused = realtype(0.0)
	hlast = realtype(0.0)
	hcur = realtype(0.0)
	tcur = realtype(0.0)
	ret = ida.IDAGetIntegratorStats(ida_mem.obj, ctypes.byref(nsteps), ctypes.byref(nrevals), ctypes.byref(nlinsetups), ctypes.byref(netfails), ctypes.byref(qlast), ctypes.byref(qcur), ctypes.byref(hinused), ctypes.byref(hlast), ctypes.byref(hcur), ctypes.byref(tcur))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetIntegratorStats() failed with flag %i"%(ret))
	return (nsteps.value, nrevals.value, nlinsetups.value, netfails.value, qlast.value, qcur.value, hinused.value, hlast.value, hcur.value, tcur.value)
ida.IDAGetIntegratorStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(realtype ), ctypes.POINTER(realtype ), ctypes.POINTER(realtype ), ctypes.POINTER(realtype )]
ida.IDAGetIntegratorStats.restype = ctypes.c_int

def IDAGetNumNonlinSolvIters(ida_mem):
	nniters = ctypes.c_long(0)
	ret = ida.IDAGetNumNonlinSolvIters(ida_mem.obj, ctypes.byref(nniters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumNonlinSolvIters() failed with flag %i"%(ret))
	return nniters.value
ida.IDAGetNumNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumNonlinSolvIters.restype = ctypes.c_int

def IDAGetNumNonlinSolvConvFails(ida_mem):
	nncfails = ctypes.c_long(0)
	ret = ida.IDAGetNumNonlinSolvConvFails(ida_mem.obj, ctypes.byref(nncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNumNonlinSolvConvFails() failed with flag %i"%(ret))
	return nncfails.value
ida.IDAGetNumNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNumNonlinSolvConvFails.restype = ctypes.c_int

def IDAGetNonlinSolvStats(ida_mem):
	nniters = ctypes.c_long(0)
	nncfails = ctypes.c_long(0)
	ret = ida.IDAGetNonlinSolvStats(ida_mem.obj, ctypes.byref(nniters), ctypes.byref(nncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDAGetNonlinSolvStats() failed with flag %i"%(ret))
	return (nniters.value, nncfils.value)
ida.IDAGetNonlinSolvStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDAGetNonlinSolvStats.restype = ctypes.c_int

def IDAGetReturnFlagName(flag):
	return ida.IDAGetReturnFlagName(flag)
ida.IDAGetReturnFlagName.argtypes = [ctypes.c_int]
ida.IDAGetReturnFlagName.restype = ctypes.c_char_p

#def IDAFree(ida_mem):
#	ret = ida.IDAFree(ida_mem.obj)
#	if ret < 0:
#		raise AssertionError("SUNDIALS ERROR: IDAFree() failed with flag %i"%(ret))

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
	if (func == None):
		return ctypes.cast(None, ATimesFn)
	exec 'def __CallbackInterface_%s(A_data, v, z):\n\treturn __ActualCallback[%i](A_data, nvecserial.NVector(v), nvecserial.NVector(z))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = ATimesFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

PSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int)
def WrapCallbackPSolveFn(func):
	if (func == None):
		return ctypes.cast(None, PSolveFn)
	exec 'def __CallbackInterface_%s(P_data, r, z, lr):\n\treturn __ActualCallback[%i](P_data, nvecserial.NVector(r), nvecserial.NVector(z), lr)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = PSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def ModifiedGS(v, h, k, p, new_vk_norm):
	ret = ida.ModifiedGS(v.data, h, k, p, new_vk_norm)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ModifiedGS() failed with flag %i"%(ret))
ida.ModifiedGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype)]
ida.ModifiedGS.restype = ctypes.c_int

def ClassicalGS(v, h, k, p, new_vk_norm, temp, s):
	ret = ida.ClassicalGS(v.data, h, k, p, new_vk_norm, temp.data, s)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ClassicalGS() failed with flag %i"%(ret))
ida.ClassicalGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype)]
ida.ClassicalGS.restype = ctypes.c_int

def QRfact(n, h, q, job):
	ret = ida.QRfact(n, h, q, job)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRfact() failed with flag %i"%(ret))
ida.QRfact.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.c_int]
ida.QRfact.restype = ctypes.c_int

def QRsol(n, h, q, b):
	ret = ida.QRsol(n, h, q, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRsol() failed with flag %i"%(ret))
ida.QRsol.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
ida.QRsol.restype = ctypes.c_int


#########################
# sundials_smalldense.h #
#########################

def denalloc(m, n):
	ret = ida.denalloc(m, n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denalloc could not allocate memory')
ida.denalloc.argtypes = [ctypes.c_long, ctypes.c_long]
ida.denalloc.restype = ctypes.POINTER(ctypes.POINTER(realtype))

def denallocpiv(n):
	ret = ida.denallocpiv(n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denallocpiv could not allocate memory')
ida.denallocpiv.argtypes = [ctypes.c_long]
ida.denallocpiv.restype = ctypes.POINTER(ctypes.c_long)

def denGETRF(a, m, n, p):
	return ida.denGETRF(a, m, n, p)
ida.denGETRF.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]
ida.denGETRF.restype = ctypes.c_long

def denGETRS(a, n, p, b):
	return ida.denGETRS(a, n, p, b)
ida.denGETRS.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
ida.denGETRS.restype = None

def denzero(a, m, n):
	ida.denzero(a, m, n)
ida.denzero.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
ida.denzero.restype = None

def dencopy(a, b, m, n):
	ida.dencopy(a, b, m, n)
ida.dencopy.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
ida.dencopy.restype = None

def denscale(c, a, m, n):
	ida.denscale(c, a, m, n)
ida.denscale.argtypes = [realtype, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
ida.denscale.restype = None

def denaddI(a, n):
	ida.denaddI(a, n)
ida.denaddI.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long]
ida.denaddI.restype = None

def denfreepiv(p):
	ida.denfreepiv(p)
ida.denfreepiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
ida.denfreepiv.restype = None

def denfree(a):
	ida.denfree(a)
ida.denfree.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype))]
ida.denfree.restype = None

def denprint(a, m, n):
	ida.denprint(a, m, n)
ida.denprint.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
ida.denprint.restype = None

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
	return ida.SpbcgMalloc(l_max, vec_tmpl.data)
ida.SpbcgMalloc.argtypes = [ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
ida.SpbcgMalloc.restype = SpbcgMem

def SpbcgSolve(mem, A_data, x, b, pretype, delta, P_data, sx, sb, atimes, psolve, res_norm, nli, nps):
	ret = ida.SpbcgSolve(mem, A_data, x.data, b.data, pretype, delta, P_data, sx.data, sb.data, atimes, psolve, res_norm, nli, nps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgSolve() failed with flag %i"%(ret))
ida.SpbcgSolve.argtypes = [SpbcgMem, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ATimesFn, PSolveFn, ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
ida.SpbcgSolve.restype = ctypes.c_int

def SpbcgFree(mem):
	ret = ida.SpbcgFree(mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgFree() failed with flag %i"%(ret))
ida.SpbcgFree.argtypes = [SpbcgMem]
ida.SpbcgFree.restype = None


##############
# ida_band.h #
##############

#IDABAND return values

IDABAND_SUCCESS = 0
IDABAND_MEM_NULL = -1
IDABAND_LMEM_NULL = -2
IDABAND_ILL_INPUT = -3
IDABAND_MEM_FAIL = -4

#Additional last_flag values

IDABAND_JACFUNC_UNRECVR = -5
IDABAND_JACFUNC_RECVR = -6

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
		"""Instantiates a new dense representation of a banded matrix.\n\tinit\tdimension of matrix to create; always square [int]\n\tmu\twidth of portion of band above diagonal\n\tml\twidth of portion of band below diagonal\n\tsmu\tis the storage upper bandwidth, mu <= smu <= size-1.  The BandGBTRF routine writes the LU factors into the storage for A. The upper triangular factor U, however, may have an upper bandwidth as big as MIN(size-1,mu+ml) because of partial pivoting. The smu field holds the upper bandwidth allocated for A."""
		if type(init) == int:
			self.size = init
			self.mu = mu
			self.ml = ml
			if smu < mu:
				raise ValueError("smu must be greater than or equal to mu")
			self.smu = smu
			self.data = ida.BandAllocMat(init, mu, ml, smu)
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
				ida.BandFreeMat(self.data)
			except:
				pass

def BandAllocMat(N, mu, ml, smu):
	"""Allocates memory for a banded matrix. Should not be called directly, rather instatiate a BandMat object."""
	return BandMat(N, mu, ml, smu)
ida.BandAllocMat.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long]
ida.BandAllocMat.restype = ctypes.POINTER(_BandMat)

def BandAllocPiv(N):
	"""BandAllocPiv allocates memory for pivot information to be filled in by the BandGBTRF routine during the factorization of an N by N band matrix. Returns a pointer, which should be passed as is to the [CV]BandPiv* functions.\n\tN\tthe size of the banded matrix for which to allocate pivot information space [int]"""
	return ida.BandAllocPiv(N)
ida.BandAllocPiv.argtypes = [ctypes.c_int]
ida.BandAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def BandGBTRF(A, p):
	"""BandGBTRF performs the LU factorization of the N by N band matrix A. This is done using standard Gaussian elimination with partial pivoting.\n\nA successful LU factorization leaves the "matrix" A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower triangular matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower triangular part of A contains the multipliers, I-L.\n\nBandGBTRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero.\n\nImportant Note: A must be allocated to accommodate the increase in upper bandwidth that occurs during factorization. If mathematically, A is a band matrix with upper bandwidth mu and lower bandwidth ml, then the upper triangular factor U can have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower triangular factor L has lower bandwidth ml. Allocate A with call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are as defined above. The user does not have to zero the "extra" storage allocated for the purpose of factorization. This will handled by the BandGBTRF routine.  ret = ida.BandGBTRF(A, p)"""
	ret = ida.BandGBTRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRF() failed with flag %i"%(ret))
ida.BandGBTRF.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long)]
ida.BandGBTRF.restype = ctypes.c_long

def BandGBTRS(A, p, b):
	"""BandGBTRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in BandGBTRF. The solution x is returned in b. This routine cannot fail if the corresponding call to BandGBTRF did not fail."""
	ret = ida.BandGBTRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRS() failed with flag %i"%(ret))
ida.BandGBTRS.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
ida.BandGBTRS.restype = None

def BandZero(A):
	"""Zeroes the entire matrix\n\tA\tthe matrix [BandMat]"""
	ida.BandZero(A.data)
ida.BandZero.argtypes = [_BandMat] 
ida.BandZero.restype = None

def BandCopy(A, B, copymu, copyml):
	"""Copies the contents of banded matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [BandMat]\n\tB\tthe source matrix [BandMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	ida.BandCopy(A.data, B.data, copymu, copyml)
ida.BandCopy.argtypes = [_BandMat, _BandMat, ctypes.c_long, ctypes.c_long]
ida.BandCopy.restype = None

def BandScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [BandMat]\n\tc\tthe scale factor [float]"""
	ida.BandScale(c, A.data)
ida.BandScale.argtypes = [realtype, _BandMat]
ida.BandScale.restype = None

def BandAddI(A):
	"""Adds 1.0 to the diagonal of the banded matrix A\n\tA\tthe matrix [BandMat]"""
	ida.BandAddI(A.data)
ida.BandAddI.argtypes = [_BandMat]
ida.BandAddI.restype = None

def BandFreeMat(A):
	"""BandFreeMat frees the memory allocated by BandAllocMat for the band matrix A.\n\tA\tthe matrix [BandMat]"""
	ida.BandFreeMat(A.data)
	del A.data
ida.BandFreeMat.argtypes = [_BandMat]
ida.BandFreeMat.restype = None

def BandFreePiv(p):
	"""BandFreePiv frees the pivot information storage memory p allocated by BandAllocPiv."""
	ida.BandFreePiv(p)
ida.BandFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
ida.BandFreePiv.restype = None

def BandPrint(A):
	"""Print out the banded matix A"""
	ida.BandPrint(A.data)
ida.BandPrint
ida.BandPrint.restype = None

IDABandJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.c_void_p, ctypes.POINTER(_BandMat), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackIDABandJacFn(func):
	if func == None:
		return ctypes.cast(None, IDABandJacFn)
	exec 'def __CallbackInterface_%s(Neq, mupper, mlower, tt, yy, yp, rr, c_j, jac_data, Jac, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](Neq, mupper, mlower, tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), c_j, jac_data, BandMat(Jac), nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDABandJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def IDABand(ida_mem, Neq, mupper, mlower):
	ret = ida.IDABand(ida_mem.obj, Neq, mupper, mlower)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABand() failed with flag %i"%(ret))
ida.IDABand.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
ida.IDABand.restype = ctypes.c_int

def IDABandSetJacFn(ida_mem, bjac, jac_data):
	ret = ida.IDABandSetJacFn(ida_mem.obj, WrapCallbackIDABandJacFn(bjac), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABandSetJacFn() failed with flag %i"%(ret))
ida.IDABandSetJacFn.argtypes = [ctypes.c_void_p, IDABandJacFn, ctypes.c_void_p]
ida.IDABandSetJacFn.restype = ctypes.c_int

def IDABandGetWorkSpace(ida_mem):
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = ida.IDABandGetWorkSpace(ida_mem.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABandGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
ida.IDABandGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDABandGetWorkSpace.restype = ctypes.c_int

def IDABandGetNumJacEvals(ida_mem):
	njevals = ctypes.c_long(0)
	ret = ida.IDABandGetNumJacEvals(ida_mem.obj, ctypes.byref(njevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABandGetNumJacEvals() failed with flag %i"%(ret))
	return njevals.value
ida.IDABandGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDABandGetNumJacEvals.restype = ctypes.c_int

def IDABandGetNumResEvals(ida_mem):
	nrevalsLS = ctypes.c_long(0)
	ret = ida.IDABandGetNumResEvals(ida_mem.obj, ctypes.byref(nrevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABandGetNumResEvals() failed with flag %i"%(ret))
	return nrevalsLS.value
ida.IDABandGetNumResEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDABandGetNumResEvals.restype = ctypes.c_int

def IDABandGetLastFlag(ida_mem):
	flag = ctypes.c_int(0)
	ret = ida.IDABandGetLastFlag(ida_mem.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABandGetLastFlag() failed with flag %i"%(ret))
	return flag.value
ida.IDABandGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDABandGetLastFlag.restype = ctypes.c_int

def IDABandGetReturnFlagName(flag):
	return ida.IDABandGetReturnFlagName(flag)
ida.IDABandGetReturnFlagName.argtypes = [ctypes.c_int]
ida.IDABandGetReturnFlagName.restype = ctypes.c_char_p


################
# ida_bbdpre.h #
################

#IDABBDPRE return values

IDABBDPRE_SUCCESS = 0
IDABBDPRE_PDATA_NULL = -11
IDABBDPRE_FUNC_UNRECVR = -12

IDABBDLocalFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackIDABBDLocalFn(func):
	if func == None:
		return ctypes.cast(None, IDABBDLocalFn)
	exec 'def __CallbackInterface_%s(Nlocal, tt, yy, yp, gval, res_data):\n\treturn __ActualCallback[%i](Nlocal, tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(gval), res_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDABBDLocalFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDABBDCommFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackIDABBDCommFn(func):
	if func == None:
		return ctypes.cast(None, IDABBDCommFn)
	exec 'def __CallbackInterface_%s(Nlocal, tt, yy, yp, res_data):\n\treturn __ActualCallback[%i](Nlocal, tt, nvecserial.NVector(yy), nvecserial.NVector(yp), res_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDABBDCommFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def IDABBDPrecAlloc(ida_mem, Nlocal, mudq, mldq, mukeep, mlkeep, dq_rel_yy, Gres, Gcomm):
	ret = ida.IDABBDPrecAlloc(ida_mem.obj, Nlocal, mudq, mldq, mukeep, mlkeep, dq_rel_yy, WrapCallbackIDABBDLocalFn(Gres), WrapCallbackIDABBDCommFn(Gcomm))
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: IDABBDPrecAlloc() failed to allocate memory")
	return ret
ida.IDABBDPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, IDABBDLocalFn, IDABBDCommFn]
ida.IDABBDPrecAlloc.restype = ctypes.c_void_p

def IDABBDSptfqmr(ida_mem, maxl, bbd_data):
	ret = ida.IDABBDSptfqmr(ida_mem.obj, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDSptfqmr() failed with flag %i"%(ret))
ida.IDABBDSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
ida.IDABBDSptfqmr.restype = ctypes.c_int

def IDABBDSpbcg(ida_mem, maxl, bbd_data):
	ret = ida.IDABBDSpbcg(ida_mem.obj, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDSpbcg() failed with flag %i"%(ret))
ida.IDABBDSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
ida.IDABBDSpbcg.restype = ctypes.c_int

def IDABBDSpgmr(ida_mem, maxl, bbd_data):
	ret = ida.IDABBDSpgmr(ida_mem.obj, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDSpgmr() failed with flag %i"%(ret))
ida.IDABBDSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
ida.IDABBDSpgmr.restype = ctypes.c_int

def IDABBDPrecReInit(bbd_data, mudq, mldq, dq_rel_yy, Gres, Gcomm):
	ret = ida.IDABBDPrecReInit(bbd_data, mudq, mldq, dq_rel_yy, WrapCallbackIDABBDLocalFn(Gres), WrapCallbackIDABBDCommFn(Gcomm))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDPrecReInit() failed with flag %i"%(ret))
ida.IDABBDPrecReInit.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, realtype, IDABBDLocalFn, IDABBDCommFn]
ida.IDABBDPrecReInit.restype = ctypes.c_int

def IDABBDPrecFree(bbd_data):
	ret = ida.IDABBDPrecFree(bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDPrecFree() failed with flag %i"%(ret))
ida.IDABBDPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
ida.IDABBDPrecFree.restype = None

def IDABBDPrecGetWorkSpace(bbd_data):
	lenrwBBDP = ctypes.c_long(0)
	leniwBBDP = ctypes.c_long(0)
	ret = ida.IDABBDPrecGetWorkSpace(bbd_data, ctypes.byref(lenrwBBDP), ctypes.byref(leniwBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwBBDP.value, leniwBBDP.value)
ida.IDABBDPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDABBDPrecGetWorkSpace.restype = ctypes.c_int

def IDABBDPrecGetNumGfnEvals(bbd_data):
	ngevalsBBDP = ctypes.c_long(0)
	ret = ida.IDABBDPrecGetNumGfnEvals(bbd_data, ctypes.byref(ngevalsBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDABBDPrecGetNumGfnEvals() failed with flag %i"%(ret))
	return ngevalsBBDP.value
ida.IDABBDPrecGetNumGfnEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDABBDPrecGetNumGfnEvals.restype = ctypes.c_int

def IDABBDPrecGetReturnFlagName(flag):
	return ida.IDABBDPrecGetReturnFlagName(flag)
ida.IDABBDPrecGetReturnFlagName.argtypes = [ctypes.c_int]
ida.IDABBDPrecGetReturnFlagName.restype = ctypes.c_char_p


###############
# ida_dense.h #
###############

#IDADENSE return values

IDADENSE_SUCCESS = 0
IDADENSE_MEM_NULL = -1
IDADENSE_LMEM_NULL = -2
IDADENSE_ILL_INPUT = -3
IDADENSE_MEM_FAIL = -4

#Additional last_flag values

IDADENSE_JACFUNC_UNRECVR = -5
IDADENSE_JACFUNC_RECVR = -6

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
		"""Instantiates a dense matrix object\n\tM\tnumber of rows\n\tN\tnumber of columns"""
		if type(M) == int and N is not None and type(N) == int:
			self.M = N
			self.N = M
			self.data = ida.DenseAllocMat(N,M)
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
				ida.DenseFreeMat(self.data)
			except:
				pass

def DenseAllocMat(M, N):
	"""Allocates memory for a dense matrix. Should not be called directly, rather instatiate a DenseMat object."""
	return DenseMat(M, N)
ida.DenseAllocMat.argtypes = [ctypes.c_int, ctypes.c_int]
ida.DenseAllocMat.restype = ctypes.POINTER(_DenseMat)

def DenseAllocPiv(N):
	"""DenseAllocPiv allocates memory for pivot information to be filled in by the DenseGETRF routine during the factorization of an N by N dense matrix. Returns a pointer, which should be passed as is to the [CV]DensePiv* functions.\n\tN\tthe size of the dense matrix for which to allocate pivot information space [int]"""
	return ida.DenseAllocPiv(N)
ida.DenseAllocPiv.argtypes = [ctypes.c_int]
ida.DenseAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def DenseGETRF(A, p):
	"""DenseGETRF performs the LU factorization of the M by N dense matrix A. This is done using standard Gaussian elimination with partial (row) pivoting. Note that this only applies to matrices with M >= N and full column rank.\n\nA successful LU factorization leaves the matrix A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower trapezoidal part of A contains the multipliers, I-L.\n\nFor square matrices (M=N), L is unit lower triangular.\n\nDenseGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	ret = ida.DenseGETRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRF() failed with flag %i"%(ret))
ida.DenseGETRF.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long)]
ida.DenseGETRF.restype = ctypes.c_long

def DenseGETRS(A, p, b):
	"""DenseGETRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in DenseGETRF. The solution x is returned in b. This routine cannot fail if the corresponding call to DenseGETRF did not fail.\nDenseGETRS does NOT check for a squre matrix!"""
	ret = ida.DenseGETRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRS() failed with flag %i"%(ret))
ida.DenseGETRS.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
ida.DenseGETRS.restype = None

def DenseZero(A):
	"""Zeroes the entire matrix\n\tA\tthe matrix [DenseMat]"""
	ida.DenseZero(A.data)
ida.DenseZero.argtypes = [_DenseMat] 
ida.DenseZero.restype = None

def DenseCopy(A, B):
	"""Copies the contents of dense matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [DenseMat]\n\tB\tthe source matrix [DenseMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	ida.DenseCopy(A.data, B.data)
ida.DenseCopy.argtypes = [_DenseMat, _DenseMat]
ida.DenseCopy.restype = None

def DenseScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [DenseMat]\n\tc\tthe scale factor [float]"""
	ida.DenseScale(c, A.data)
ida.DenseScale.argtypes = [realtype, _DenseMat]
ida.DenseScale.restype = None

def DenseAddI(A):
	"""Adds 1.0 to the diagonal of the dense matrix A\n\tA\tthe matrix [DenseMat]"""
	ida.DenseAddI(A.data)
ida.DenseAddI.argtypes = [_DenseMat]
ida.DenseAddI.restype = None

def DenseFreeMat(A):
	"""DenseFreeMat frees the memory allocated by DenseAllocMat for the band matrix A.\n\tA\tthe matrix [DenseMat]"""
	ida.DenseFreeMat(A.data)
	del A.data
ida.DenseFreeMat.argtypes = [_DenseMat]
ida.DenseFreeMat.restype = None

def DenseFreePiv(p):
	"""DenseFreePiv frees the pivot information storage memory p allocated by DenseAllocPiv."""
	ida.DenseFreePiv(p)
ida.DenseFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
ida.DenseFreePiv.restype = None

def DensePrint(A):
	"""Print out the dense matix A"""
	ida.DensePrint(A.data)
ida.DensePrint
ida.DensePrint.restype = None

IDADenseJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.c_void_p, ctypes.POINTER(_DenseMat), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackIDADenseJacFn(func):
	if func == None:
		return ctypes.cast(None, IDADenseJacFn)
	exec 'def __CallbackInterface_%s(Neq, tt, yy, yp, rr, c_j, jac_data, Jac, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](Neq, tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), c_j, jac_data, DenseMat(Jac), nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDADenseJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def IDADense(ida_mem, Neq):
	ret = ida.IDADense(ida_mem.obj, Neq)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADense() failed with flag %i"%(ret))
ida.IDADense.argtypes = [ctypes.c_void_p, ctypes.c_long]
ida.IDADense.restype = ctypes.c_int

def IDADenseSetJacFn(ida_mem, djac, jac_data):
	ret = ida.IDADenseSetJacFn(ida_mem.obj, WrapCallbackIDADenseJacFn(djac), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADenseSetJacFn() failed with flag %i"%(ret))
ida.IDADenseSetJacFn.argtypes = [ctypes.c_void_p, IDADenseJacFn, ctypes.c_void_p]
ida.IDADenseSetJacFn.restype = ctypes.c_int

def IDADenseGetWorkSpace(ida_mem):
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = ida.IDADenseGetWorkSpace(ida_mem.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADenseGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
ida.IDADenseGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDADenseGetWorkSpace.restype = ctypes.c_int

def IDADenseGetNumJacEvals(ida_mem):
	njevals = ctypes.c_long(0)
	ret = ida.IDADenseGetNumJacEvals(ida_mem.obj, ctypes.byref(njevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADenseGetNumJacEvals() failed with flag %i"%(ret))
	return njevals.value
ida.IDADenseGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDADenseGetNumJacEvals.restype = ctypes.c_int

def IDADenseGetNumResEvals(ida_mem):
	nrevalsLS = ctypes.c_long(0)
	ret = ida.IDADenseGetNumResEvals(ida_mem.obj, ctypes.byref(nrevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADenseGetNumResEvals() failed with flag %i"%(ret))
	return nrevalsLS.value
ida.IDADenseGetNumResEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDADenseGetNumResEvals.restype = ctypes.c_int

def IDADenseGetLastFlag(ida_mem):
	flag = ctypes.c_int(0)
	ret = ida.IDADenseGetLastFlag(ida_mem.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDADenseGetLastFlag() failed with flag %i"%(ret))
	return flag.value
ida.IDADenseGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDADenseGetLastFlag.restype = ctypes.c_int

def IDADenseGetReturnFlagName(flag):
	return ida.IDADenseGetReturnFlagName(flag)
ida.IDADenseGetReturnFlagName.argtypes = [ctypes.c_int]
ida.IDADenseGetReturnFlagName.restype = ctypes.c_char_p


###############
# ida_spbcg.h #
###############

def IDASpbcg(ida_mem, maxl):
	ret = ida.IDASpbcg(ida_mem.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpbcg() failed with flag %i"%(ret))
ida.IDASpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASpbcg.restype = ctypes.c_int


###############
# ida_spgmr.h #
###############

def IDASpgmr(ida_mem, maxl):
	ret = ida.IDASpgmr(ida_mem.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpgmr() failed with flag %i"%(ret))
ida.IDASpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASpgmr.restype = ctypes.c_int


###############
# ida_spils.h #
###############

#IDASPILS return values 

IDASPILS_SUCCESS = 0
IDASPILS_MEM_NULL = -1
IDASPILS_LMEM_NULL = -2
IDASPILS_ILL_INPUT = -3
IDASPILS_MEM_FAIL = -4

IDASpilsPrecSetupFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackIDASpilsPrecSetupFn(func):
	if func == None:
		return ctypes.cast(None, IDASpilsPrecSetupFn)
	exec 'def __CallbackInterface_%s(tt, yy, yp, rr, c_j, prec_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), c_j, prec_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDASpilsPrecSetupFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDASpilsPrecSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackIDASpilsPrecSolveFn(func):
	if func == None:
		return ctypes.cast(None, IDASpilsPrecSolveFn)
	exec 'def __CallbackInterface_%s(tt, yy, yp, rr, rvec, zvec, c_j, delta, prec_data, tmp):\n\treturn __ActualCallback[%i](tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), nvecserial.NVector(rvec), nvecserial.NVector(zvec), c_j, delta, prec_data, nvecserial.NVector(tmp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDASpilsPrecSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

IDASpilsJacTimesVecFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackIDASpilsJacTimesVecFn(func):
	if func == None:
		return ctypes.cast(None, IDASpilsJacTimesVecFn)
	exec 'def __CallbackInterface_%s(tt, yy, yp, rr, v, Jv, c_j, jac_data, tmp1, tmp2):\n\treturn __ActualCallback[%i](tt, nvecserial.NVector(yy), nvecserial.NVector(yp), nvecserial.NVector(rr), nvecserial.NVector(v), nvecserial.NVector(Jv), c_j, jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = IDASpilsJacTimesVecFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def IDASpilsSetPreconditioner(ida_mem, pset, psolve, prec_data):
	ret = ida.IDASpilsSetPreconditioner(ida_mem.obj, WrapCallbackIDASpilsPrecSetupFn(pset), WrapCallbackIDASpilsPrecSolveFn(psolve), prec_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetPreconditioner() failed with flag %i"%(ret))
ida.IDASpilsSetPreconditioner.argtypes = [ctypes.c_void_p, IDASpilsPrecSetupFn, IDASpilsPrecSolveFn, ctypes.c_void_p]
ida.IDASpilsSetPreconditioner.restype = ctypes.c_int

def IDASpilsSetJacTimesVecFn(ida_mem, jtimes, jac_data):
	ret = ida.IDASpilsSetJacTimesVecFn(ida_mem.obj, WrapCallbackIDASpilsJacTimesVecFn(jtimes), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetJacTimesVecFn() failed with flag %i"%(ret))
ida.IDASpilsSetJacTimesVecFn.argtypes = [ctypes.c_void_p, IDASpilsJacTimesVecFn, ctypes.c_void_p]
ida.IDASpilsSetJacTimesVecFn.restype = ctypes.c_int

def IDASpilsSetGSType(ida_mem, gstype):
	ret = ida.IDASpilsSetGSType(ida_mem.obj, gstype)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetGSType() failed with flag %i"%(ret))
ida.IDASpilsSetGSType.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASpilsSetGSType.restype = ctypes.c_int

def IDASpilsSetMaxRestarts(ida_mem, maxrs):
	ret = ida.IDASpilsSetMaxRestarts(ida_mem.obj, maxrs)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetMaxRestarts() failed with flag %i"%(ret))
ida.IDASpilsSetMaxRestarts.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASpilsSetMaxRestarts.restype = ctypes.c_int

def IDASpilsSetMaxl(ida_mem, maxl):
	ret = ida.IDASpilsSetMaxl(ida_mem.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetMaxl() failed with flag %i"%(ret))
ida.IDASpilsSetMaxl.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASpilsSetMaxl.restype = ctypes.c_int

def IDASpilsSetEpsLin(ida_mem, eplifac):
	ret = ida.IDASpilsSetEpsLin(ida_mem.obj, eplifac)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetEpsLin() failed with flag %i"%(ret))
ida.IDASpilsSetEpsLin.argtypes = [ctypes.c_void_p, realtype]
ida.IDASpilsSetEpsLin.restype = ctypes.c_int

def IDASpilsSetIncrementFactor(ida_mem, dqincfac):
	ret = ida.IDASpilsSetIncrementFactor(ida_mem.obj, dqincfac)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsSetIncrementFactor() failed with flag %i"%(ret))
ida.IDASpilsSetIncrementFactor.argtypes = [ctypes.c_void_p, realtype]
ida.IDASpilsSetIncrementFactor.restype = ctypes.c_int

def IDASpilsGetWorkSpace(ida_mem):
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = ida.IDASpilsGetWorkSpace(ida_mem.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
ida.IDASpilsGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetWorkSpace.restype = ctypes.c_int

def IDASpilsGetNumPrecEvals(ida_mem):
	npevals = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumPrecEvals(ida_mem.obj, ctypes.byref(npevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumPrecEvals() failed with flag %i"%(ret))
	return npevals.value
ida.IDASpilsGetNumPrecEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumPrecEvals.restype = ctypes.c_int

def IDASpilsGetNumPrecSolves(ida_mem):
	npsolves = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumPrecSolves(ida_mem.obj, ctypes.byref(npsolves))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumPrecSolves() failed with flag %i"%(ret))
	return npsolves.value
ida.IDASpilsGetNumPrecSolves.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumPrecSolves.restype = ctypes.c_int

def IDASpilsGetNumLinIters(ida_mem):
	nliters = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumLinIters(ida_mem.obj, ctypes.byref(nliters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumLinIters() failed with flag %i"%(ret))
	return nliters.value
ida.IDASpilsGetNumLinIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumLinIters.restype = ctypes.c_int

def IDASpilsGetNumConvFails(ida_mem):
	nlcfails = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumConvFails(ida_mem.obj, ctypes.byref(nlcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumConvFails() failed with flag %i"%(ret))
	return nlcfails.value
ida.IDASpilsGetNumConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumConvFails.restype = ctypes.c_int

def IDASpilsGetNumJtimesEvals(ida_mem):
	njvevals = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumJtimesEvals(ida_mem.obj, ctypes.byref(njvevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumJtimesEvals() failed with flag %i"%(ret))
	return njvevals.value
ida.IDASpilsGetNumJtimesEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumJtimesEvals.restype = ctypes.c_int

def IDASpilsGetNumResEvals(ida_mem):
	nrevalsLS = ctypes.c_long(0)
	ret = ida.IDASpilsGetNumResEvals(ida_mem.obj, ctypes.byref(nrevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetNumResEvals() failed with flag %i"%(ret))
	return nrevalsLS.value
ida.IDASpilsGetNumResEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
ida.IDASpilsGetNumResEvals.restype = ctypes.c_int

def IDASpilsGetLastFlag(ida_mem):
	flag = ctypes.c_int(0)
	ret = ida.IDASpilsGetLastFlag(ida_mem.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASpilsGetLastFlag() failed with flag %i"%(ret))
	return flag.value
ida.IDASpilsGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
ida.IDASpilsGetLastFlag.restype = ctypes.c_int

def IDASpilsGetReturnFlagName(flag):
	return ida.IDASpilsGetReturnFlagName(flag)
ida.IDASpilsGetReturnFlagName.argtypes = [ctypes.c_int]
ida.IDASpilsGetReturnFlagName.restype = ctypes.c_char_p


#################
# ida_sptfqmr.h #
#################

def IDASptfqmr(ida_mem, maxl):
	ret = ida.IDASptfqmr(ida_mem.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: IDASptfqmr() failed with flag %i"%(ret))
ida.IDASptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int]
ida.IDASptfqmr.restype = ctypes.c_int
