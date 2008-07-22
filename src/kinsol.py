"""Python bindings for the kinsol, kinsol_band, kinsol_bbdpre, kinsol_dense, kinsol_spbcgs, kinsol_spmgr, kinsol_spils, and kinsol_sptfqmr header files."""
import ctypes
import sundials_core
import nvecserial

realtype = nvecserial.realtype
NVector = nvecserial.NVector

kinsol = sundials_core.loadlib("kinsol")

#KINSOL return flags
KIN_SUCCESS = 0
KIN_INITIAL_GUESS_OK = 1
KIN_STEP_LT_STPTOL = 2

KIN_WARNING = 99

KIN_MEM_NULL = -1
KIN_ILL_INPUT = -2
KIN_NO_MALLOC = -3
KIN_MEM_FAIL = -4
KIN_LINESEARCH_NONCONV = -5
KIN_MAXITER_REACHED = -6
KIN_MXNEWT_5X_EXCEEDED = -7
KIN_LINESEARCH_BCFAIL = -8
KIN_LINSOLV_NO_RECOVERY = -9
KIN_LINIT_FAIL = -10
KIN_LSETUP_FAIL = -11
KIN_LSOLVE_FAIL = -12

KIN_SYSFUNC_FAIL = -13
KIN_FIRST_SYSFUNC_ERR = -14
KIN_REPTD_SYSFUNC_ERR = -15

#Enumeration for inputs to KINSetEtaForm
KIN_ETACHOICE1 = 1
KIN_ETACHOICE2 = 2
KIN_ETACONSTANT = 3

#Enumeration for global strategy
KIN_NONE = 0
KIN_LINESEARCH = 1

__Callback = []
__ActualCallback = []

class KinsolMemObj(object):
	def __init__(self, obj):
		self.obj = obj

	def __del__(self):
		kinsol.KINFree(self.obj)

KINSysFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackKINSysFn(func):
	if func == None:
		return ctypes.cast(None, KINSysFn)
	exec 'def __CallbackInterface_%s(uu, fval, f_data):\n\treturn __ActualCallback[%i](nvecserial.NVector(uu), nvecserial.NVector(fval), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINSysFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

KINErrHandlerFn = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_void_p)
def WrapCallbackKINErrHandlerFn(func):
	if func == None:
		return ctypes.cast(None, KINErrHandlerFn)
	exec 'def __CallbackInterface_%s(error_code, module, function, msg, eh_data):\n\treturn __ActualCallback[%i](error_code, module, function, msg, eh_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINErrHandlerFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

KINInfoHandlerFn = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_void_p)
def WrapCallbackKINInfoHandlerFn(func):
	if func == None:
		return ctypes.cast(None, KINInfoHandlerFn)
	exec 'def __CallbackInterface_%s(module, function, msg, ih_data):\n\treturn __ActualCallback[%i](module, function, msg, ih_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINInfoHandlerFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def KINCreate():
	obj = kinsol.KINCreate()
	if obj == None:
		raise AssertionError("SUNDIALS ERROR: KINCreate() failed - returned NULL pointer")
	return KinsolMemObj(obj)
kinsol.KINCreate.argtypes = []
kinsol.KINCreate.restype = ctypes.c_void_p

def KINSetErrHandlerFn(kinsolmemobj, ehfun, eh_data):
	ret = kinsol.KINSetErrHandlerFn(kinsolmemobj.obj, WrapCallbackKINErrHandlerFn(ehfun), eh_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetErrHandlerFn() failed with flag %i"%(ret))
kinsol.KINSetErrHandlerFn.argtypes = [ctypes.c_void_p, KINErrHandlerFn, ctypes.c_void_p]
kinsol.KINSetErrHandlerFn.restype = ctypes.c_int

def KINSetErrFile(kinsolmemobj, errfp):
	ret = kinsol.KINSetErrFile(kinsolmemobj.obj, sundials_core.fdopen(errfp.fileno, errfp.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetErrFile() failed with flag %i"%(ret))
kinsol.KINSetErrFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
kinsol.KINSetErrFile.restype = ctypes.c_int

def KINSetInfoHandlerFn(kinsolmemobj, ihfun, ih_data):
	ret = kinsol.KINSetInfoHandlerFn(kinsolmemobj.obj, WrapCallbackKINInfoHandlerFn(ihfun), ih_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetInfoHandlerFn() failed with flag %i"%(ret))
kinsol.KINSetInfoHandlerFn.argtypes = [ctypes.c_void_p, KINInfoHandlerFn, ctypes.c_void_p]
kinsol.KINSetInfoHandlerFn.restype = ctypes.c_int

def KINSetInfoFile(kinsolmemobj, infofp):
	ret = kinsol.KINSetInfoFile(kinsolmemobj.obj, sundials_core.fdopen(infofp.fileno, infofp.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetInfoFile() failed with flag %i"%(ret))
kinsol.KINSetInfoFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
kinsol.KINSetInfoFile.restype = ctypes.c_int

def KINSetFdata(kinsolmemobj, f_data):
	ret = kinsol.KINSetFdata(kinsolmemobj.obj, f_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetFdata() failed with flag %i"%(ret))
kinsol.KINSetFdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
kinsol.KINSetFdata.restype = ctypes.c_int

def KINSetPrintLevel(kinmemm, printfl):
	ret = kinsol.KINSetPrintLevel(kinmemm, printfl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetPrintLevel() failed with flag %i"%(ret))
kinsol.KINSetPrintLevel.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSetPrintLevel.restype = ctypes.c_int

def KINSetNumMaxIters(kinsolmemobj, mxiter):
	ret = kinsol.KINSetNumMaxIters(kinsolmemobj.obj, mxiter)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetNumMaxIters() failed with flag %i"%(ret))
kinsol.KINSetNumMaxIters.argtypes = [ctypes.c_void_p, ctypes.c_long]
kinsol.KINSetNumMaxIters.restype = ctypes.c_int

def KINSetNoInitSetup(kinsolmemobj, noInitSetup):
	ret = kinsol.KINSetNoInitSetup(kinsolmemobj.obj, ctypes.byref(noInitSetup))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetNoInitSetup() failed with flag %i"%(ret))
kinsol.KINSetNoInitSetup.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSetNoInitSetup.restype = ctypes.c_int

def KINSetNoResMon(kinsolmemobj, noNNIResMon):
	ret = kinsol.KINSetNoResMon(kinsolmemobj.obj, ctypes.byref(noNNIResMon))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetNoResMon() failed with flag %i"%(ret))
kinsol.KINSetNoResMon.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSetNoResMon.restype = ctypes.c_int

def KINSetMaxSetupCalls(kinsolmemobj, msbset):
	ret = kinsol.KINSetMaxSetupCalls(kinsolmemobj.obj, msbset)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetMaxSetupCalls() failed with flag %i"%(ret))
kinsol.KINSetMaxSetupCalls.argtypes = [ctypes.c_void_p, ctypes.c_long]
kinsol.KINSetMaxSetupCalls.restype = ctypes.c_int

def KINSetMaxSubSetupCalls(kinsolmemobj, msbsetsub):
	ret = kinsol.KINSetMaxSubSetupCalls(kinsolmemobj.obj, msbsetsub)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetMaxSubSetupCalls() failed with flag %i"%(ret))
kinsol.KINSetMaxSubSetupCalls.argtypes = [ctypes.c_void_p, ctypes.c_long]
kinsol.KINSetMaxSubSetupCalls.restype = ctypes.c_int

def KINSetEtaForm(kinsolmemobj, etachoice):
	ret = kinsol.KINSetEtaForm(kinsolmemobj.obj, etachoice)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetEtaForm() failed with flag %i"%(ret))
kinsol.KINSetEtaForm.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSetEtaForm.restype = ctypes.c_int

def KINSetEtaConstValue(kinsolmemobj, eta):
	ret = kinsol.KINSetEtaConstValue(kinsolmemobj.obj, eta)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetEtaConstValue() failed with flag %i"%(ret))
kinsol.KINSetEtaConstValue.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetEtaConstValue.restype = ctypes.c_int

def KINSetEtaParams(kinsolmemobj, egamma, ealpha):
	ret = kinsol.KINSetEtaParams(kinsolmemobj.obj, egamma, ealpha)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetEtaParams() failed with flag %i"%(ret))
kinsol.KINSetEtaParams.argtypes = [ctypes.c_void_p, realtype, realtype]
kinsol.KINSetEtaParams.restype = ctypes.c_int

def KINSetResMonParams(kinsolmemobj, omegamin, omegamax):
	ret = kinsol.KINSetResMonParams(kinsolmemobj.obj, omegamin, omegamax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetResMonParams() failed with flag %i"%(ret))
kinsol.KINSetResMonParams.argtypes = [ctypes.c_void_p, realtype, realtype]
kinsol.KINSetResMonParams.restype = ctypes.c_int

def KINSetResMonConstValue(kinsolmemobj, omegaconst):
	ret = kinsol.KINSetResMonConstValue(kinsolmemobj.obj, omegaconst)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetResMonConstValue() failed with flag %i"%(ret))
kinsol.KINSetResMonConstValue.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetResMonConstValue.restype = ctypes.c_int

def KINSetNoMinEps(kinsolmemobj, noMinEps):
	ret = kinsol.KINSetNoMinEps(kinsolmemobj.obj, ctypes.byref(noMinEps))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetNoMinEps() failed with flag %i"%(ret))
kinsol.KINSetNoMinEps.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSetNoMinEps.restype = ctypes.c_int

def KINSetMaxNewtonStep(kinsolmemobj, mxnewtstep):
	ret = kinsol.KINSetMaxNewtonStep(kinsolmemobj.obj, mxnewtstep)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetMaxNewtonStep() failed with flag %i"%(ret))
kinsol.KINSetMaxNewtonStep.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetMaxNewtonStep.restype = ctypes.c_int

def KINSetMaxBetaFails(kinsolmemobj, mxnbcf):
	ret = kinsol.KINSetMaxBetaFails(kinsolmemobj.obj, mxnbcf)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetMaxBetaFails() failed with flag %i"%(ret))
kinsol.KINSetMaxBetaFails.argtypes = [ctypes.c_void_p, ctypes.c_long]
kinsol.KINSetMaxBetaFails.restype = ctypes.c_int

def KINSetRelErrFunc(kinsolmemobj, relfunc):
	ret = kinsol.KINSetRelErrFunc(kinsolmemobj.obj, relfunc)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetRelErrFunc() failed with flag %i"%(ret))
kinsol.KINSetRelErrFunc.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetRelErrFunc.restype = ctypes.c_int

def KINSetFuncNormTol(kinsolmemobj, fnormtol):
	ret = kinsol.KINSetFuncNormTol(kinsolmemobj.obj, fnormtol)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetFuncNormTol() failed with flag %i"%(ret))
kinsol.KINSetFuncNormTol.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetFuncNormTol.restype = ctypes.c_int

def KINSetScaledStepTol(kinsolmemobj, scsteptol):
	ret = kinsol.KINSetScaledStepTol(kinsolmemobj.obj, scsteptol)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetScaledStepTol() failed with flag %i"%(ret))
kinsol.KINSetScaledStepTol.argtypes = [ctypes.c_void_p, realtype]
kinsol.KINSetScaledStepTol.restype = ctypes.c_int

def KINSetConstraints(kinsolmemobj, constraints):
	ret = kinsol.KINSetConstraints(kinsolmemobj.obj, constraints.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetConstraints() failed with flag %i"%(ret))
kinsol.KINSetConstraints.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
kinsol.KINSetConstraints.restype = ctypes.c_int

def KINSetSysFunc(kinsolmemobj, func):
	ret = kinsol.KINSetSysFunc(kinsolmemobj.obj, WrapCallbackKINSysFn(func))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSetSysFunc() failed with flag %i"%(ret))
kinsol.KINSetSysFunc.argtypes = [ctypes.c_void_p, KINSysFn]
kinsol.KINSetSysFunc.restype = ctypes.c_int

def KINMalloc(kinsolmemobj, func, tmpl):
	ret = kinsol.KINMalloc(kinsolmemobj.obj, WrapCallbackKINSysFn(func), tmpl.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINMalloc() failed with flag %i"%(ret))
kinsol.KINMalloc.argtypes = [ctypes.c_void_p, KINSysFn, ctypes.POINTER(nvecserial._NVector)]
kinsol.KINMalloc.restype = ctypes.c_int

def KINSol(kinsolmemobj, uu, strategy, u_scale, f_scale):
	ret = kinsol.KINSol(kinsolmemobj.obj, uu.data, strategy, u_scale.data, f_scale.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSol() failed with flag %i"%(ret))
kinsol.KINSol.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector)]
kinsol.KINSol.restype = ctypes.c_int

def KINGetWorkSpace(kinsolmemobj):
	lenrw = ctypes.c_long(0)
	leniw = ctypes.c_long(0)
	ret = kinsol.KINGetWorkSpace(kinsolmemobj.obj, ctypes.byref(lenrw), ctypes.byref(leniw))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetWorkSpace() failed with flag %i"%(ret))
	return (lenrw.value, leniw.value)
kinsol.KINGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
kinsol.KINGetWorkSpace.restype = ctypes.c_int

def KINGetNumNonlinSolvIters(kinsolmemobj):
	nniters = ctypes.c_long(0)
	ret = kinsol.KINGetNumNonlinSolvIters(kinsolmemobj.obj, ctypes.byref(nniters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetNumNonlinSolvIters() failed with flag %i"%(ret))
	return nniters.value
kinsol.KINGetNumNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINGetNumNonlinSolvIters.restype = ctypes.c_int

def KINGetNumFuncEvals(kinsolmemobj):
	nfevals = ctypes.c_long(0)
	ret = kinsol.KINGetNumFuncEvals(kinsolmemobj.obj, ctypes.byref(nfevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetNumFuncEvals() failed with flag %i"%(ret))
	return nfevals.value
kinsol.KINGetNumFuncEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINGetNumFuncEvals.restype = ctypes.c_int

def KINGetNumBetaCondFails(kinsolmemobj):
	nbcfails = ctypes.c_long(0)
	ret = kinsol.KINGetNumBetaCondFails(kinsolmemobj.obj, ctypes.byref(nbcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetNumBetaCondFails() failed with flag %i"%(ret))
	return nbcfails.value
kinsol.KINGetNumBetaCondFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINGetNumBetaCondFails.restype = ctypes.c_int

def KINGetNumBacktrackOps(kinsolmemobj):
	nbacktr = ctypes.c_long(0)
	ret = kinsol.KINGetNumBacktrackOps(kinsolmemobj.obj, ctypes.byref(nbacktr))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetNumBacktrackOps() failed with flag %i"%(ret))
	return nbacktr.value
kinsol.KINGetNumBacktrackOps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINGetNumBacktrackOps.restype = ctypes.c_int

def KINGetFuncNorm(kinsolmemobj):
	fnorm = realtype(0)
	ret = kinsol.KINGetFuncNorm(kinsolmemobj.obj, ctypes.byref(fnorm))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetFuncNorm() failed with flag %i"%(ret))
	return fnorm.value
kinsol.KINGetFuncNorm.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
kinsol.KINGetFuncNorm.restype = ctypes.c_int

def KINGetStepLength(kinsolmemobj):
	steplength = ctypes.c_long(0)
	ret = kinsol.KINGetStepLength(kinsolmemobj.obj, ctypes.byref(steplength))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINGetStepLength() failed with flag %i"%(ret))
	return steplength.value
kinsol.KINGetStepLength.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
kinsol.KINGetStepLength.restype = ctypes.c_int

def KINGetReturnFlagName(flag):
	return kinsol.KINGetReturnFlagName(flag)
kinsol.KINGetReturnFlagName.argtypes = [ctypes.c_int]
kinsol.KINGetReturnFlagName.restype = ctypes.c_char_p

def KINFree(kinsolmemobj):
	kinsol.KINFree(kinsolmemobj.obj)
	del kinsolmemobj.obj

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
	ret = kinsol.ModifiedGS(v.data, h, k, p, new_vk_norm)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ModifiedGS() failed with flag %i"%(ret))
kinsol.ModifiedGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype)]
kinsol.ModifiedGS.restype = ctypes.c_int

def ClassicalGS(v, h, k, p, new_vk_norm, temp, s):
	ret = kinsol.ClassicalGS(v.data, h, k, p, new_vk_norm, temp.data, s)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ClassicalGS() failed with flag %i"%(ret))
kinsol.ClassicalGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype)]
kinsol.ClassicalGS.restype = ctypes.c_int

def QRfact(n, h, q, job):
	ret = kinsol.QRfact(n, h, q, job)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRfact() failed with flag %i"%(ret))
kinsol.QRfact.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.c_int]
kinsol.QRfact.restype = ctypes.c_int

def QRsol(n, h, q, b):
	ret = kinsol.QRsol(n, h, q, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRsol() failed with flag %i"%(ret))
kinsol.QRsol.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
kinsol.QRsol.restype = ctypes.c_int


#########################
# sundials_smalldense.h #
#########################

def denalloc(m, n):
	ret = kinsol.denalloc(m, n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denalloc could not allocate memory')
kinsol.denalloc.argtypes = [ctypes.c_long, ctypes.c_long]
kinsol.denalloc.restype = ctypes.POINTER(ctypes.POINTER(realtype))

def denallocpiv(n):
	ret = kinsol.denallocpiv(n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denallocpiv could not allocate memory')
kinsol.denallocpiv.argtypes = [ctypes.c_long]
kinsol.denallocpiv.restype = ctypes.POINTER(ctypes.c_long)

def denGETRF(a, m, n, p):
	return kinsol.denGETRF(a, m, n, p)
kinsol.denGETRF.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]
kinsol.denGETRF.restype = ctypes.c_long

def denGETRS(a, n, p, b):
	return kinsol.denGETRS(a, n, p, b)
kinsol.denGETRS.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
kinsol.denGETRS.restype = None

def denzero(a, m, n):
	kinsol.denzero(a, m, n)
kinsol.denzero.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
kinsol.denzero.restype = None

def dencopy(a, b, m, n):
	kinsol.dencopy(a, b, m, n)
kinsol.dencopy.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
kinsol.dencopy.restype = None

def denscale(c, a, m, n):
	kinsol.denscale(c, a, m, n)
kinsol.denscale.argtypes = [realtype, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
kinsol.denscale.restype = None

def denaddI(a, n):
	kinsol.denaddI(a, n)
kinsol.denaddI.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long]
kinsol.denaddI.restype = None

def denfreepiv(p):
	kinsol.denfreepiv(p)
kinsol.denfreepiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
kinsol.denfreepiv.restype = None

def denfree(a):
	kinsol.denfree(a)
kinsol.denfree.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype))]
kinsol.denfree.restype = None

def denprint(a, m, n):
	kinsol.denprint(a, m, n)
kinsol.denprint.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
kinsol.denprint.restype = None

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
	return kinsol.SpbcgMalloc(l_max, vec_tmpl.data)
kinsol.SpbcgMalloc.argtypes = [ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
kinsol.SpbcgMalloc.restype = SpbcgMem

def SpbcgSolve(mem, A_data, x, b, pretype, delta, P_data, sx, sb, atimes, psolve, res_norm, nli, nps):
	ret = kinsol.SpbcgSolve(mem, A_data, x.data, b.data, pretype, delta, P_data, sx.data, sb.data, atimes, psolve, res_norm, nli, nps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgSolve() failed with flag %i"%(ret))
kinsol.SpbcgSolve.argtypes = [SpbcgMem, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ATimesFn, PSolveFn, ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
kinsol.SpbcgSolve.restype = ctypes.c_int

def SpbcgFree(mem):
	ret = kinsol.SpbcgFree(mem)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgFree() failed with flag %i"%(ret))
kinsol.SpbcgFree.argtypes = [SpbcgMem]
kinsol.SpbcgFree.restype = None


#################
# kinsol_band.h #
#################

#KINBAND return values

KINBAND_SUCCESS = 0
KINBAND_MEM_NULL = -1
KINBAND_LMEM_NULL = -2
KINBAND_ILL_INPUT = -3
KINBAND_MEM_FAIL = -4

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
			self.data = kinsol.BandAllocMat(init, mu, ml, smu)
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
				kinsol.BandFreeMat(self.data)
			except:
				pass

def BandAllocMat(N, mu, ml, smu):
	"""Allocates memory for a banded matrix. Should not be called directly, rather instatiate a BandMat object."""
	return BandMat(N, mu, ml, smu)
kinsol.BandAllocMat.argtypes = [ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long]
kinsol.BandAllocMat.restype = ctypes.POINTER(_BandMat)

def BandAllocPiv(N):
	"""BandAllocPiv allocates memory for pivot information to be filled in by the BandGBTRF routine during the factorization of an N by N band matrix. Returns a pointer, which should be passed as is to the [CV]BandPiv* functions.\n\tN\tthe size of the banded matrix for which to allocate pivot information space [int]"""
	return kinsol.BandAllocPiv(N)
kinsol.BandAllocPiv.argtypes = [ctypes.c_int]
kinsol.BandAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def BandGBTRF(A, p):
	"""BandGBTRF performs the LU factorization of the N by N band matrix A. This is done using standard Gaussian elimination with partial pivoting.\n\nA successful LU factorization leaves the "matrix" A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower triangular matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower triangular part of A contains the multipliers, I-L.\n\nBandGBTRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero.\n\nImportant Note: A must be allocated to accommodate the increase in upper bandwidth that occurs during factorization. If mathematically, A is a band matrix with upper bandwidth mu and lower bandwidth ml, then the upper triangular factor U can have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower triangular factor L has lower bandwidth ml. Allocate A with call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are as defined above. The user does not have to zero the "extra" storage allocated for the purpose of factorization. This will handled by the BandGBTRF routine.  ret = kinsol.BandGBTRF(A, p)"""
	ret = kinsol.BandGBTRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRF() failed with flag %i"%(ret))
kinsol.BandGBTRF.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long)]
kinsol.BandGBTRF.restype = ctypes.c_long

def BandGBTRS(A, p, b):
	"""BandGBTRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in BandGBTRF. The solution x is returned in b. This routine cannot fail if the corresponding call to BandGBTRF did not fail."""
	ret = kinsol.BandGBTRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: BandGBTRS() failed with flag %i"%(ret))
kinsol.BandGBTRS.argtypes = [_BandMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
kinsol.BandGBTRS.restype = None

def BandZero(A):
	"""Zeroes the entire matrix\n\tA\tthe matrix [BandMat]"""
	kinsol.BandZero(A.data)
kinsol.BandZero.argtypes = [_BandMat] 
kinsol.BandZero.restype = None

def BandCopy(A, B, copymu, copyml):
	"""Copies the contents of banded matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [BandMat]\n\tB\tthe source matrix [BandMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	kinsol.BandCopy(A.data, B.data, copymu, copyml)
kinsol.BandCopy.argtypes = [_BandMat, _BandMat, ctypes.c_long, ctypes.c_long]
kinsol.BandCopy.restype = None

def BandScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [BandMat]\n\tc\tthe scale factor [float]"""
	kinsol.BandScale(c, A.data)
kinsol.BandScale.argtypes = [realtype, _BandMat]
kinsol.BandScale.restype = None

def BandAddI(A):
	"""Adds 1.0 to the diagonal of the banded matrix A\n\tA\tthe matrix [BandMat]"""
	kinsol.BandAddI(A.data)
kinsol.BandAddI.argtypes = [_BandMat]
kinsol.BandAddI.restype = None

def BandFreeMat(A):
	"""BandFreeMat frees the memory allocated by BandAllocMat for the band matrix A.\n\tA\tthe matrix [BandMat]"""
	kinsol.BandFreeMat(A.data)
	del A.data
kinsol.BandFreeMat.argtypes = [_BandMat]
kinsol.BandFreeMat.restype = None

def BandFreePiv(p):
	"""BandFreePiv frees the pivot information storage memory p allocated by BandAllocPiv."""
	kinsol.BandFreePiv(p)
kinsol.BandFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
kinsol.BandFreePiv.restype = None

def BandPrint(A):
	"""Print out the banded matix A"""
	kinsol.BandPrint(A.data)
kinsol.BandPrint
kinsol.BandPrint.restype = None

KINBandJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.POINTER(_BandMat), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackKINBandJacFn(func):
	if func == None:
		return ctypes.cast(None, KINBandJacFn)
	exec 'def __CallbackInterface_%s(N, mupper, mlower, J, u, fu, jac_data, tmp1, tmp2):\n\treturn __ActualCallback[%i](N, mupper, mlower, BandMat(J), nvecserial.NVector(u), nvecserial.NVector(fu), jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINBandJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def KINBand(kinsolmemobj, N, mupper, mlower):
	ret = kinsol.KINBand(kinsolmemobj.obj, N, mupper, mlower)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBand() failed with flag %i"%(ret))
kinsol.KINBand.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
kinsol.KINBand.restype = ctypes.c_int

def KINBandSetJacFn(kinsolmemobj, bjac, jac_data):
	ret = kinsol.KINBandSetJacFn(kinsolmemobj.obj, WrapCallbackKINBandJacFn(bjac), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBandSetJacFn() failed with flag %i"%(ret))
kinsol.KINBandSetJacFn.argtypes = [ctypes.c_void_p, KINBandJacFn, ctypes.c_void_p]
kinsol.KINBandSetJacFn.restype = ctypes.c_int

def KINBandGetWorkSpace(kinsolmemobj):
	lenrwB = ctypes.c_long(0)
	leniwB = ctypes.c_long(0)
	ret = kinsol.KINBandGetWorkSpace(kinsolmemobj.obj, ctypes.byref(lenrwB), ctypes.byref(leniwB))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBandGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwB.value, leniwB.value)
kinsol.KINBandGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
kinsol.KINBandGetWorkSpace.restype = ctypes.c_int

def KINBandGetNumJacEvals(kinsolmemobj):
	njevalsB = ctypes.c_long(0)
	ret = kinsol.KINBandGetNumJacEvals(kinsolmemobj.obj, ctypes.byref(njevalsB))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBandGetNumJacEvals() failed with flag %i"%(ret))
	return njevalsB.value
kinsol.KINBandGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINBandGetNumJacEvals.restype = ctypes.c_int

def KINBandGetNumFuncEvals(kinsolmemobj):
	nfevalsB = ctypes.c_long(0)
	ret = kinsol.KINBandGetNumFuncEvals(kinsolmemobj.obj, ctypes.byref(nfevalsB))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBandGetNumFuncEvals() failed with flag %i"%(ret))
	return nfevalsB.value
kinsol.KINBandGetNumFuncEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINBandGetNumFuncEvals.restype = ctypes.c_int

def KINBandGetLastFlag(kinsolmemobj):
	flag = ctypes.c_long(0)
	ret = kinsol.KINBandGetLastFlag(kinsolmemobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBandGetLastFlag() failed with flag %i"%(ret))
	return flag.value
kinsol.KINBandGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
kinsol.KINBandGetLastFlag.restype = ctypes.c_int

def KINBandGetReturnFlagName(flag):
	return kinsol.KINBandGetReturnFlagName(flag)
kinsol.KINBandGetReturnFlagName.argtypes = [ctypes.c_int]
kinsol.KINBandGetReturnFlagName.restype = ctypes.c_char_p


###################
# kinsol_bbdpre.h #
###################

#KINBBDPRE return values

KINBBDPRE_SUCCESS = 0
KINBBDPRE_PDATA_NULL = -11
KINBBDPRE_FUNC_UNRECVR = -12

KINCommFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackKINCommFn(func):
	if func == None:
		return ctypes.cast(None, KINCommFn)
	exec 'def __CallbackInterface_%s(Nlocal, u, f_data):\n\treturn __ActualCallback[%i](Nlocal, nvecserial.NVector(u), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINCommFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

KINLocalFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackKINLocalFn(func):
	if func == None:
		return ctypes.cast(None, KINLocalFn)
	exec 'def __CallbackInterface_%s(Nlocal, uu, gval, f_data):\n\treturn __ActualCallback[%i](Nlocal, nvecserial.NVector(uu), nvecserial.NVector(gval), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINLocalFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def KINBBDPrecAlloc(kinsolmemobj, Nlocal, mudq, mldq, mukeep, mlkeep, dq_rel_uu, gloc, gcomm):
	ret = kinsol.KINBBDPrecAlloc(kinsolmemobj.obj, Nlocal, mudq, mldq, mukeep, mlkeep, dq_rel_uu, WrapCallbackKINLocalFn(gloc), WrapCallbackKINCommFn(gcomm))
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecAlloc() failed to allocate memory")
	return ret
kinsol.KINBBDPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, KINLocalFn, KINCommFn]
kinsol.KINBBDPrecAlloc.restype = ctypes.c_void_p

def KINBBDSptfqmr(kinsolmemobj, maxl, p_data):
	ret = kinsol.KINBBDSptfqmr(kinsolmemobj.obj, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDSptfqmr() failed with flag %i"%(ret))
kinsol.KINBBDSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
kinsol.KINBBDSptfqmr.restype = ctypes.c_int

def KINBBDSpbcg(kinsolmemobj, maxl, p_data):
	ret = kinsol.KINBBDSpbcg(kinsolmemobj.obj, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDSpbcg() failed with flag %i"%(ret))
kinsol.KINBBDSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
kinsol.KINBBDSpbcg.restype = ctypes.c_int

def KINBBDSpgmr(kinsolmemobj, maxl, p_data):
	ret = kinsol.KINBBDSpgmr(kinsolmemobj.obj, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDSpgmr() failed with flag %i"%(ret))
kinsol.KINBBDSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
kinsol.KINBBDSpgmr.restype = ctypes.c_int

def KINBBDPrecFree(p_data):
	ret = kinsol.KINBBDPrecFree(p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecFree() failed with flag %i"%(ret))
kinsol.KINBBDPrecFree.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
kinsol.KINBBDPrecFree.restype = None

def KINBBDPrecGetWorkSpace(p_data):
	lenrwBBDP = ctypes.c_long(0)
	leniwBBDP = ctypes.c_long(0)
	ret = kinsol.KINBBDPrecGetWorkSpace(p_data, ctypes.byref(lenrwBBDP), ctypes.byref(leniwBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwBBDP.value, leniwBBDP.value)
kinsol.KINBBDPrecGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
kinsol.KINBBDPrecGetWorkSpace.restype = ctypes.c_int

def KINBBDPrecGetNumGfnEvals(p_data):
	ngevalsBBDP = ctypes.c_long(0)
	ret = kinsol.KINBBDPrecGetNumGfnEvals(p_data, ctypes.byref(ngevalsBBDP))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecGetNumGfnEvals() failed with flag %i"%(ret))
	return ngevalsBBDP.value
kinsol.KINBBDPrecGetNumGfnEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINBBDPrecGetNumGfnEvals.restype = ctypes.c_int

def KINBBDPrecGetReturnFlagName(flag):
	return kinsol.KINBBDPrecGetReturnFlagName(flag)
kinsol.KINBBDPrecGetReturnFlagName.argtypes = [ctypes.c_int]
kinsol.KINBBDPrecGetReturnFlagName.restype = ctypes.c_char_p

def KINBBDPrecSetup(uu, uscale, fval, fscale, p_data, vtemp1, vtemp2):
	ret = kinsol.KINBBDPrecSetup(uu.data, uscale.data, fval.data, fscale.data, p_data, vtemp1.data, vtemp2.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecSetup() failed with flag %i"%(ret))
kinsol.KINBBDPrecSetup.argtypes = [ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector)]
kinsol.KINBBDPrecSetup.restype = ctypes.c_int

def KINBBDPrecSolve(uu, uscale, fval, fscale, vv, p_data, vtemp):
	ret = kinsol.KINBBDPrecSolve(uu.data, uscale.data, fval.data, fscale.data, vv.data, p_data, vtemp.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINBBDPrecSolve() failed with flag %i"%(ret))
kinsol.KINBBDPrecSolve.argtypes = [ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
kinsol.KINBBDPrecSolve.restype = ctypes.c_int


##################
# kinsol_dense.h #
##################

#KINDENSE return values

KINDENSE_SUCCESS = 0
KINDENSE_MEM_NULL = -1
KINDENSE_LMEM_NULL = -2
KINDENSE_ILL_INPUT = -3
KINDENSE_MEM_FAIL = -4

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
			self.data = kinsol.DenseAllocMat(N,M)
			self.copy = False
		elif type(M) == ctypes.POINTER(_DenseMat):
			self.M = M.contents.M
			self.N = M.contents.N
			self.data = M
			self.copy = True
		else:
			raise TypeError("Cannot initialize DenseMat from type %s"%(type(init).__name__))
	
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
				kinsol.DenseFreeMat(self.data)
			except:
				pass

def DenseAllocMat(M, N):
	"""Allocates memory for a dense matrix. Should not be called directly, rather instatiate a DenseMat object."""
	return DenseMat(M, N)
kinsol.DenseAllocMat.argtypes = [ctypes.c_int, ctypes.c_int]
kinsol.DenseAllocMat.restype = ctypes.POINTER(_DenseMat)

def DenseAllocPiv(N):
	"""DenseAllocPiv allocates memory for pivot information to be filled in by the DenseGETRF routine during the factorization of an N by N dense matrix. Returns a pointer, which should be passed as is to the [CV]DensePiv* functions.\n\tN\tthe size of the dense matrix for which to allocate pivot information space [int]"""
	return kinsol.DenseAllocPiv(N)
kinsol.DenseAllocPiv.argtypes = [ctypes.c_int]
kinsol.DenseAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def DenseGETRF(A, p):
	"""DenseGETRF performs the LU factorization of the M by N dense matrix A. This is done using standard Gaussian elimination with partial (row) pivoting. Note that this only applies to matrices with M >= N and full column rank.\n\nA successful LU factorization leaves the matrix A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower trapezoidal part of A contains the multipliers, I-L.\n\nFor square matrices (M=N), L is unit lower triangular.\n\nDenseGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	ret = kinsol.DenseGETRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRF() failed with flag %i"%(ret))
kinsol.DenseGETRF.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long)]
kinsol.DenseGETRF.restype = ctypes.c_long

def DenseGETRS(A, p, b):
	"""DenseGETRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in DenseGETRF. The solution x is returned in b. This routine cannot fail if the corresponding call to DenseGETRF did not fail.\nDenseGETRS does NOT check for a squre matrix!"""
	ret = kinsol.DenseGETRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRS() failed with flag %i"%(ret))
kinsol.DenseGETRS.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
kinsol.DenseGETRS.restype = None

def DenseZero(A):
	"""Zeroes the entire matrix\n\tA\tthe matrix [DenseMat]"""
	kinsol.DenseZero(A.data)
kinsol.DenseZero.argtypes = [_DenseMat] 
kinsol.DenseZero.restype = None

def DenseCopy(A, B):
	"""Copies the contents of dense matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [DenseMat]\n\tB\tthe source matrix [DenseMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	kinsol.DenseCopy(A.data, B.data)
kinsol.DenseCopy.argtypes = [_DenseMat, _DenseMat]
kinsol.DenseCopy.restype = None

def DenseScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [DenseMat]\n\tc\tthe scale factor [float]"""
	kinsol.DenseScale(c, A.data)
kinsol.DenseScale.argtypes = [realtype, _DenseMat]
kinsol.DenseScale.restype = None

def DenseAddI(A):
	"""Adds 1.0 to the diagonal of the dense matrix A\n\tA\tthe matrix [DenseMat]"""
	kinsol.DenseAddI(A.data)
kinsol.DenseAddI.argtypes = [_DenseMat]
kinsol.DenseAddI.restype = None

def DenseFreeMat(A):
	"""DenseFreeMat frees the memory allocated by DenseAllocMat for the band matrix A.\n\tA\tthe matrix [DenseMat]"""
	kinsol.DenseFreeMat(A.data)
	del A.data
kinsol.DenseFreeMat.argtypes = [_DenseMat]
kinsol.DenseFreeMat.restype = None

def DenseFreePiv(p):
	"""DenseFreePiv frees the pivot information storage memory p allocated by DenseAllocPiv."""
	kinsol.DenseFreePiv(p)
kinsol.DenseFreePiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
kinsol.DenseFreePiv.restype = None

def DensePrint(A):
	"""Print out the dense matix A"""
	kinsol.DensePrint(A.data)
kinsol.DensePrint
kinsol.DensePrint.restype = None

KINDenseJacFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, ctypes.POINTER(_DenseMat), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackKINDenseJacFn(func):
	if func == None:
		return ctypes.cast(None, KINDenseJacFn)
	exec 'def __CallbackInterface_%s(N, J, uu, fval, jac_data, vtemp1, vtemp2):\n\treturn __ActualCallback[%i](N, DenseMat(J), nvecserial.NVector(uu), nvecserial.NVector(fval), jac_data, nvecserial.NVector(vtemp1), nvecserial.NVector(vtemp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINDenseJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def KINDense(kinsolmemobj, N):
	ret = kinsol.KINDense(kinsolmemobj.obj, N)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDense() failed with flag %i"%(ret))
kinsol.KINDense.argtypes = [ctypes.c_void_p, ctypes.c_long]
kinsol.KINDense.restype = ctypes.c_int

def KINDenseSetJacFn(kinsolmemobj, djac, jac_data):
	ret = kinsol.KINDenseSetJacFn(kinsolmemobj.obj, WrapCallbackKINDenseJacFn(djac), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDenseSetJacFn() failed with flag %i"%(ret))
kinsol.KINDenseSetJacFn.argtypes = [ctypes.c_void_p, KINDenseJacFn, ctypes.c_void_p]
kinsol.KINDenseSetJacFn.restype = ctypes.c_int

def KINDenseGetWorkSpace(kinsolmemobj):
	lenrwD = ctypes.c_long(0)
	leniwD = ctypes.c_long(0)
	ret = kinsol.KINDenseGetWorkSpace(kinsolmemobj.obj, ctypes.byref(lenrwD), ctypes.byref(leniwD))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDenseGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwD.value, leniwD.value)
kinsol.KINDenseGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
kinsol.KINDenseGetWorkSpace.restype = ctypes.c_int

def KINDenseGetNumJacEvals(kinsolmemobj):
	njevalsD = ctypes.c_long(0)
	ret = kinsol.KINDenseGetNumJacEvals(kinsolmemobj.obj, ctypes.byref(njevalsD))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDenseGetNumJacEvals() failed with flag %i"%(ret))
	return njevalsD.value
kinsol.KINDenseGetNumJacEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINDenseGetNumJacEvals.restype = ctypes.c_int

def KINDenseGetNumFuncEvals(kinsolmemobj):
	nfevalsD = ctypes.c_long(0)
	ret = kinsol.KINDenseGetNumFuncEvals(kinsolmemobj.obj, ctypes.byref(nfevalsD))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDenseGetNumFuncEvals() failed with flag %i"%(ret))
	return nfevalsD.value
kinsol.KINDenseGetNumFuncEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINDenseGetNumFuncEvals.restype = ctypes.c_int

def KINDenseGetLastFlag(kinsolmemobj):
	flag = ctypes.c_int(0)
	ret = kinsol.KINDenseGetLastFlag(kinsolmemobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINDenseGetLastFlag() failed with flag %i"%(ret))
	return flag.value
kinsol.KINDenseGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
kinsol.KINDenseGetLastFlag.restype = ctypes.c_int

def KINDenseGetReturnFlagName(flag):
	return kinsol.KINDenseGetReturnFlagName(flag)
kinsol.KINDenseGetReturnFlagName.argtypes = [ctypes.c_int]
kinsol.KINDenseGetReturnFlagName.restype = ctypes.c_char_p


##################
# kinsol_spbcg.h #
##################

def KINSpbcg(kinsolmemobj, maxl):
	ret = kinsol.KINSpbcg(kinsolmemobj.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpbcg() failed with flag %i"%(ret))
kinsol.KINSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSpbcg.restype = ctypes.c_int


##################
# kinsol_spgmr.h #
##################

def KINSpgmr(kinsolmemobj, maxl):
	ret = kinsol.KINSpgmr(kinsolmemobj.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpgmr() failed with flag %i"%(ret))
kinsol.KINSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSpgmr.restype = ctypes.c_int


##################
# kinsol_spils.h #
##################

#KINSPILS return values

KINSPILS_SUCCESS = 0

KINSPILS_MEM_NULL = -1
KINSPILS_LMEM_NULL = -2
KINSPILS_ILL_INPUT = -3
KINSPILS_MEM_FAIL = -4

# KINSPILS solver constant

KINSPILS_MAXL = 10

KINSpilsPrecSetupFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector))
def WrapCallbackKINSpilsPrecSetupFn(func):
	if func == None:
		return ctypes.cast(None, KINSpilsPrecSetupFn)
	exec 'def __CallbackInterface_%s(uu, uscale, fval, fscale, P_data, vtemp1, vtemp2):\n\treturn __ActualCallback[%i](nvecserial.NVector(uu), nvecserial.NVector(uscale), nvecserial.NVector(fval), nvecserial.NVector(fscale), P_data, nvecserial.NVector(vtemp1), nvecserial.NVector(vtemp2))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINSpilsPrecSetupFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

KINSpilsPrecSolveFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector))
def WrapCallbackKINSpilsPrecSolveFn(func):
	if func == None:
		return ctypes.cast(None, KINSpilsPrecSolveFn)
	exec 'def __CallbackInterface_%s(uu, uscale, fval, fscale, vv, P_data, vtemp):\n\treturn __ActualCallback[%i](nvecserial.NVector(uu), nvecserial.NVector(uscale), nvecserial.NVector(fval), nvecserial.NVector(fscale), nvecserial.NVector(vv), P_data, nvecserial.NVector(vtemp))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINSpilsPrecSolveFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

KINSpilsJacTimesVecFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
def WrapCallbackKINSpilsJacTimesVecFn(func):
	if func == None:
		return ctypes.cast(None, KINSpilsJacTimesVecFn)
	exec 'def __CallbackInterface_%s(v, Jv, uu, new_uu, J_data):\n\treturn __ActualCallback[%i](nvecserial.NVector(v), nvecserial.NVector(Jv), nvecserial.NVector(uu), ctypes.byref(new_uu), J_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = KINSpilsJacTimesVecFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def KINSpilsSetMaxRestarts(kinsolmemobj, maxrs):
	ret = kinsol.KINSpilsSetMaxRestarts(kinsolmemobj.obj, maxrs)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsSetMaxRestarts() failed with flag %i"%(ret))
kinsol.KINSpilsSetMaxRestarts.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSpilsSetMaxRestarts.restype = ctypes.c_int

def KINSpilsSetPreconditioner(kinsolmemobj, pset, psolve, P_data):
	ret = kinsol.KINSpilsSetPreconditioner(kinsolmemobj.obj, WrapCallbackKINSpilsPrecSetupFn(pset), WrapCallbackKINSpilsPrecSolveFn(psolve), P_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsSetPreconditioner() failed with flag %i"%(ret))
kinsol.KINSpilsSetPreconditioner.argtypes = [ctypes.c_void_p, KINSpilsPrecSetupFn, KINSpilsPrecSolveFn, ctypes.c_void_p]
kinsol.KINSpilsSetPreconditioner.restype = ctypes.c_int

def KINSpilsSetJacTimesVecFn(kinsolmemobj, jtimes, J_data):
	ret = kinsol.KINSpilsSetJacTimesVecFn(kinsolmemobj.obj, WrapCallbackKINSpilsJacTimesVecFn(jtimes), J_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsSetJacTimesVecFn() failed with flag %i"%(ret))
kinsol.KINSpilsSetJacTimesVecFn.argtypes = [ctypes.c_void_p, KINSpilsJacTimesVecFn, ctypes.c_void_p]
kinsol.KINSpilsSetJacTimesVecFn.restype = ctypes.c_int

def KINSpilsGetWorkSpace(kinsolmemobj):
	lenrwSG = ctypes.c_long(0)
	leniwSG = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetWorkSpace(kinsolmemobj.obj, ctypes.byref(lenrwSG), ctypes.byref(leniwSG))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwSG.value, leniwSG.value)
kinsol.KINSpilsGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetWorkSpace.restype = ctypes.c_int

def KINSpilsGetNumPrecEvals(kinsolmemobj):
	npevals = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumPrecEvals(kinsolmemobj.obj, ctypes.byref(npevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumPrecEvals() failed with flag %i"%(ret))
	return npevals.value
kinsol.KINSpilsGetNumPrecEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumPrecEvals.restype = ctypes.c_int

def KINSpilsGetNumPrecSolves(kinsolmemobj):
	npsolves = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumPrecSolves(kinsolmemobj.obj, ctypes.byref(npsolves))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumPrecSolves() failed with flag %i"%(ret))
	return npsolves.value
kinsol.KINSpilsGetNumPrecSolves.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumPrecSolves.restype = ctypes.c_int

def KINSpilsGetNumLinIters(kinsolmemobj):
	nliters = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumLinIters(kinsolmemobj.obj, ctypes.byref(nliters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumLinIters() failed with flag %i"%(ret))
	return nliters.value
kinsol.KINSpilsGetNumLinIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumLinIters.restype = ctypes.c_int

def KINSpilsGetNumConvFails(kinsolmemobj):
	nlcfails = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumConvFails(kinsolmemobj.obj, ctypes.byref(nlcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumConvFails() failed with flag %i"%(ret))
	return nlcfails.value
kinsol.KINSpilsGetNumConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumConvFails.restype = ctypes.c_int

def KINSpilsGetNumJtimesEvals(kinsolmemobj):
	njvevals = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumJtimesEvals(kinsolmemobj.obj, ctypes.byref(njvevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumJtimesEvals() failed with flag %i"%(ret))
	return njvevals.value
kinsol.KINSpilsGetNumJtimesEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumJtimesEvals.restype = ctypes.c_int

def KINSpilsGetNumFuncEvals(kinsolmemobj):
	nfevalsS = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetNumFuncEvals(kinsolmemobj.obj, ctypes.byref(nfevalsS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetNumFuncEvals() failed with flag %i"%(ret))
	return nfevalsS.value
kinsol.KINSpilsGetNumFuncEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
kinsol.KINSpilsGetNumFuncEvals.restype = ctypes.c_int

def KINSpilsGetLastFlag(kinsolmemobj):
	flag = ctypes.c_long(0)
	ret = kinsol.KINSpilsGetLastFlag(kinsolmemobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSpilsGetLastFlag() failed with flag %i"%(ret))
	return flag.value
kinsol.KINSpilsGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
kinsol.KINSpilsGetLastFlag.restype = ctypes.c_int

def KINSpilsGetReturnFlagName(flag):
	return kinsol.KINSpilsGetReturnFlagName(flag)
kinsol.KINSpilsGetReturnFlagName.argtypes = [ctypes.c_int]
kinsol.KINSpilsGetReturnFlagName.restype = ctypes.c_char_p


####################
# kinsol_sptfqmr.h #
####################

def KINSptfqmr(kinsolmemobj, maxl):
	ret = kinsol.KINSptfqmr(kinsolmemobj.obj, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: KINSptfqmr() failed with flag %i"%(ret))
kinsol.KINSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int]
kinsol.KINSptfqmr.restype = ctypes.c_int

