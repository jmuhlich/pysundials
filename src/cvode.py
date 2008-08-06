"""Python bindings for the cvode, cvode_band, cvode_bandpre, cvode_bbdpre, cvode_dense, cvode_diag, cvode_spbcgs, cvode_spmgr, cvode_spils, and cvode_sptfqmr header files."""
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
	"""Creates a wrapper around a python callable object, that can be used as a callback for the RHS function. RHS functions take time_step (float), y (NVector), ydot (NVector), and f_data (c_void_p) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVRhsFn)
	exec 'def __CallbackInterface_%s(t, y, ydot, f_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), nvecserial.NVector(ydot), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVRhsFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVRootFn = ctypes.CFUNCTYPE(ctypes.c_int, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype), ctypes.c_void_p)
def WrapCallbackCVRootFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the root finding function. Root finding functions take time_step (float), y (NVector), gout (NVector), and g_data (c_void_p) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVRootFn)
	exec 'def __CallbackInterface_%s(t, y, gout, g_data):\n\treturn __ActualCallback[%i](t, nvecserial.NVector(y), gout, g_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVRootFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVEwtFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVEwtFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the error weight function. Error weight functions take y (NVector), ewt (NVector), and e_data (c_void_p) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVEwtFn)
	exec 'def __CallbackInterface_%s(y, ewt, e_data):\n\treturn __ActualCallback[%i](nvecserial.NVector(y), nvecserial.NVector(ewt), e_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVEwtFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVErrHandlerFn = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_void_p)
def WrapCallbackCVErrHandlerFn(func):
	"""Creates a wrapper around a python callable object, that can be used as a callback for the error handler function. Error handler functions take error_code (int), module (string), function_name (string), message (string), and eh_data (c_void_p) as parameters, and have no return value."""
	if func == None:
		return ctypes.cast(None, CVErrHandlerFn)
	exec 'def __CallbackInterface_%s(error_code, module, function, msg, eh_data):\n\treturn __ActualCallback[%i](error_code, module, function, msg, eh_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVErrHandlerFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVodeCreate(lmm, iter):
	"""Creates an internal memory block for a problem to be solved by CVODE.\n\tlmm\tis the type of linear multistep method to be used. The legal values are CV_ADAMS and CV_BDF.\n\t\tCV_ADAMS\t(Adams-Moulton) is recommended for non-stiff problems.\n\t\tCV_BDF\t\t(Backward Differentiation Formula) is recommended for stiff problems.  \n\titer\tspecifies whether functional or newton iteration will be used. Legal values are: \n\t\tCV_FUNCTIONAL\tdoes not require linear algebra \n\t\tCV_NEWTON\trequires the solution of linear systems. Requires the specification of a CVODE linear solver. (Recommended for stiff problems)."""
	obj = cvode.CVodeCreate(lmm, iter)
	if obj == None:
		raise AssertionError("SUNDIALS ERROR: CVodeCreate() failed - returned NULL pointer")
	return CVodeMemObj(obj)
cvode.CVodeCreate.argtypes = [ctypes.c_int, ctypes.c_int]
cvode.CVodeCreate.restype = ctypes.c_void_p

def CVodeSetErrHandlerFn(cvodememobj, func,  eh_data):
	"""Sets a user provided error handling function.\n\tfunc\tis a python callable wrapped by WrapCallbackCVErrHandlerFn\n\teh_data\tis any user data you would like to pass"""
	ret = cvode.CVodeSetErrHandlerFn(cvodememobj.obj, WrapCallbackCVErrHandlerFn(func), eh_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrHanderFn() failed with flag %i"%(ret))
cvode.CVodeSetErrHandlerFn.argtypes = [ctypes.c_void_p, CVErrHandlerFn, ctypes.c_void_p]
cvode.CVodeSetErrHandlerFn.restype = ctypes.c_int

def CVodeSetErrFile(cvodememobj, fileobj):
	"""Sets the file where all error messages and warnings will be written if the default error handling function is used.\n\tfileobj\tis a file object opened for writing"""
	ret = cvode.CVodeSetErrFile(cvodememobj.obj, sundials_core.fdopen(fileobj.fileno, fileobj.mode))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetErrFile() failed with flag %i"%(ret))
cvode.CVodeSetErrFile.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvode.CVodeSetErrFile.restype = ctypes.c_int

def CVodeSetFdata(cvodememobj, f_data):
	"""Sets the pointer to user data. If called, whenever the user's f function is called, f_data will point to this data"""
	ret = cvode.CVodeSetFdata(cvodememobj.obj, f_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetFdata() failed with flag %i"%(ret))
cvode.CVodeSetFdata.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
cvode.CVodeSetFdata.restype = ctypes.c_int

def CVodeSetEwtFn(cvodememobj, func, e_data):
	"""Sets a user provided error weight function.\n\tfunc\tis a python callable wrapped by WrapCallbackCVEwtFn\n\te_data\tis any user data you would like to pass"""
	ret = cvode.CVodeSetEwtFn(cvodememobj.obj, WrapCallbackCVEwtFn(func), e_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetEwtFn() failed with flag %i"%(ret))
cvode.CVodeSetEwtFn.argtypes = [ctypes.c_void_p, CVEwtFn, ctypes.c_void_p]
cvode.CVodeSetEwtFn.restype = ctypes.c_int

def CVodeSetMaxOrd(cvodememobj, maxord):
	"""Sets the maximum lmm order to be used by the solver. Defaults to 12 for CV_ADAMS, and 5 for CV_BDF"""
	ret = cvode.CVodeSetMaxOrd(cvodememobj.obj, maxord)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxOrd() failed with flag %i"%(ret))
cvode.CVodeSetMaxOrd.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxOrd.restype = ctypes.c_int

def CVodeSetMaxNumSteps(cvodememobj, maxsteps):
	"""Sets the maximum number of internal steps to be taken by the solver in its attempt to reach tout."""
	ret = cvode.CVodeSetMaxNumSteps(cvodememobj.obj, maxsteps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNumSteps() failed with flag %i"%(ret))
cvode.CVodeSetMaxNumSteps.argtypes = [ctypes.c_void_p, ctypes.c_long] 
cvode.CVodeSetMaxNumSteps.restype = ctypes.c_int

def CVodeSetMaxHnilWarns(cvodememobj, mxhnil):
	"""Sets the maximum number of warning messages issued by the solver that t+h==t on the next internal step. A value of -1 means no such warnings are issued. Default is 10."""
	ret = cvode.CVodeSetMaxHnilWarns(cvodememobj.obj, mxhnil)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxHnilWarns() failed with flag %i"%(ret))
cvode.CVodeSetMaxHnilWarns.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxHnilWarns.restype = ctypes.c_int

def CVodeSetStabLimDet(cvodememobj, stldet):
	"""Turns stability limit detection on or off (stldet == 1 [ON], stldet == 0 [OFF)\nWhen CV_BDF is used and order is 3 or greater, CVsldet is called to detect stability limit. If limit is detected, the order is reduced. Default is off"""
	ret = cvode.CVodeSetStabLimDet(cvodememobj.obj, stldet)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStabLimDet() failed with flag %i"%(ret))
cvode.CVodeSetStabLimDet.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetStabLimDet.restype = ctypes.c_int

def CVodeSetInitStep(cvodememobj, hin):
	"""Sets initial step size. CVODE estimates this value by default."""
	ret = cvode.CVodeSetInitStep(cvodememobj.obj, hin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetInitStep() failed with flag %i"%(ret))
cvode.CVodeSetInitStep.argtypes = [ctypes.c_void_p, realtype] 
cvode.CVodeSetInitStep.restype = ctypes.c_int

def CVodeSetMinStep(cvodememobj, hmin):
	"""Sets the absolute minimum of step sie allowed. Default is 0.0"""
	ret = cvode.CVodeSetMinStep(cvodememobj.obj, hmin)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMinStep() failed with flag %i"%(ret))
cvode.CVodeSetMinStep.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetMinStep.restype = ctypes.c_int

def CVodeSetMaxStep(cvodememobj, hmax):
	"""Sets the maximum absolute value of step size allowed. Default is infinity."""
	ret = cvode.CVodeSetMaxStep(cvodememobj.obj, hmax)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxStep() failed with flag %i"%(ret))
cvode.CVodeSetMaxStep.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetMaxStep.restype = ctypes.c_int

def CVodeSetStopTime(cvodememobj, tstop):
	"""Sets the independent variable value past which the solution is not to proceed. Default is infinity."""
	ret = cvode.CVodeSetStopTime(cvodememobj.obj, tstop)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetStopTime() failed with flag %i"%(ret))
cvode.CVodeSetStopTime.argtypes = [ctypes.c_void_p, realtype] 
cvode.CVodeSetStopTime.restype = ctypes.c_int

def CVodeSetMaxErrTestFails(cvodememobj, maxnef):
	"""Sets the maximum number of error test failures in attempting one step. Default is 7."""
	ret = cvode.CVodeSetMaxErrTestFails(cvodememobj.obj, maxnef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxErrTestFails() failed with flag %i"%(ret))
cvode.CVodeSetMaxErrTestFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxErrTestFails.restype = ctypes.c_int

def CVodeSetMaxNonlinIters(cvodememobj, maxcor):
	"""Sets the maximum number of nonlinear solver iterations at one solution. Default is 3."""
	ret = cvode.CVodeSetMaxNonlinIters(cvodememobj.obj, maxcor)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxNonlinIters() failed with flag %i"%(ret))
cvode.CVodeSetMaxNonlinIters.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxNonlinIters.restype = ctypes.c_int

def CVodeSetMaxConvFails(cvodememobj, maxncf):
	"""Sets the maximum number of convergence failures allowed in attempting one step. Default is 10."""
	ret = cvode.CVodeSetMaxConvFails(cvodememobj.obj, maxncf)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetMaxConvFails() failed with flag %i"%(ret))
cvode.CVodeSetMaxConvFails.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetMaxConvFails.restype = ctypes.c_int

def CVodeSetNonlinConvCoef(cvodememobj, nlscoef):
	"""Sets the coefficient in the nonlinear convergence test. Default is 0.1"""
	ret = cvode.CVodeSetNonlinConvCoef(cvodememobj.obj, nlscoef)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetNonlinConvCoef() failed with flag %i"%(ret))
cvode.CVodeSetNonlinConvCoef.argtypes = [ctypes.c_void_p, realtype]
cvode.CVodeSetNonlinConvCoef.restype = ctypes.c_int

def CVodeSetIterType(cvodememobj, iter):
	"""Changes the current nonlinear iteration type. Initially set by CVodecreate()."""
	ret = cvode.CVodeSetIterType(cvodememobj.obj, iter)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetIterType() failed with flag %i"%(ret))
cvode.CVodeSetIterType.argtypes = [ctypes.c_void_p, ctypes.c_int]
cvode.CVodeSetIterType.restype = ctypes.c_int

def CVodeSetTolerances(cvodememobj, itol, reltol, abstol):
	"""Changes the integration tolerances between calls to CVode(). Initially (re)set by CVodeMalloc()/CVodeReInit()."""
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
		raise ValueError("itol must be one of CV_SS or CV_SV")
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeSetTolerances() failed with flag %i"%(ret))
cvode.CVodeSetTolerances.argtypes = [ctypes.c_void_p, ctypes.c_int, realtype, ctypes.c_void_p]
cvode.CVodeSetTolerances.restype = ctypes.c_int

def CVodeMalloc(cvodememobj, func, t0, y0, itol, reltol, abstol):
	"""CVodeMalloc allocates and initializes memory for a problem to be solved by CVODE.\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()\n\tfunc\t\ta python callable defining the side function in y' = f(t,y), as wrapped by WrapCallbackCVRhsFn\n\tt0\t\tthe initial value of t. [float]\n\ty0\t\tthe initial condition vector y(t0). [NVector]\n\titol\t\tthe type of tolerances to be used. Legal values are:\n\t\tCV_SS\tscalar relative and absolute tolerances\n\t\tCV_SV\tscalar relative tolerance and vector absolute tolerance.\n\t\tCV_WF\tindicates that the user will provide a function to evaluate the error weights. In this case, reltol and abstol are ignored.\n\treltol\t\tthe relative tolerance scalar [float]\n\tabstol\t\tabsolute tolerance scalar [float] or an N_Vector of absolute tolerances.\n\nThe parameters itol, reltol, and abstol define a vector of error weights, ewt, with components\n\tewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = CV_SS), or\n\tewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = CV_SV).\nThis vector is used in all error and convergence tests, which use a weighted RMS norm on all error-like vectors v:\n\tWRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),\nwhere N is the problem dimension."""
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
	"""CVodeReInit re-initializes CVode for the solution of a problem, where a prior call to CVodeMalloc has been made with the same problem size N. CVodeReInit performs the same input checking and initializations that CVodeMalloc does.  But it does no memory allocation, assuming that the existing internal memory is sufficient for the new problem.\n\nThe use of CVodeReInit requires that the maximum method order, maxord, is no larger for the new problem than for the problem specified in the last call to CVodeMalloc.  This condition is automatically fulfilled if the multistep method parameter lmm is unchanged (or changed from CV_ADAMS to CV_BDF) and the default value for maxord is specified."""
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
	"""CVodeRootInit initializes a rootfinding problem to be solved during the integration of the ODE system. It must be called after CVodeCreate, and before CVode. The arguments are:\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()\n\tnrtfn\t\tnumber of function whose roots are to be found [int >= 0]\n\tfunc\t\ta python callable as returned by WrapCallbackCVRootFn() which defines the functions g_i whose roots are sought.\n\tg_data\t\tuser data that will be passed to the user's g function every time g is called.\n\nIf a new problem is to be solved with a call to CVodeReInit, where the new problem has no root functions but the prior one did, then call CVodeRootInit with nrtfn = 0."""
	ret = cvode.CVodeRootInit(cvodememobj.obj, nrtfn, WrapCallbackCVRootFn(func), g_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeRootInit() failed with flag %i"%(ret))
cvode.CVodeRootInit.argtypes = [ctypes.c_void_p, ctypes.c_int, CVRootFn, ctypes.c_void_p]
cvode.CVodeRootInit.restype = ctypes.c_int

def CVode(cvodememobj, tout, yout, tret, itask):
	"""CVode integrates the ODE over an interval in t.  If itask is CV_NORMAL, then the solver integrates from its current internal t value to a point at or beyond tout, then interpolates to t = tout and returns y(tout) in the user- allocated vector yout. If itask is CV_ONE_STEP, then the solver takes one internal time step and returns in yout the value of y at the new internal time. In this case, tout is used only during the first call to CVode to determine the direction of integration and the rough scale of the t variable.  If itask is CV_NORMAL_TSTOP or CV_ONE_STEP_TSTOP, then CVode returns the solution at tstop if that comes sooner than tout or the end of the next internal step, respectively.  In any case, the time reached by the solver is placed in (*tret). The user is responsible for allocating the memory for this value.\n\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()\n\ttout\t\tthe next time at which a computed solution is desired. [float]\n\tyout\t\tthe computed solution vector. In CV_NORMAL mode with no errors and no roots found, yout=y(tout). [NVector]\n\ttret\t\tis set to the time reached by the solver [realtype]\n\titask\t\tis one of CV_NORMAL, CV_ONE_STEP, CV_NORMAL_TSTOP, or CV_ONE_STEP_TSTOP. See CVODE documentation for more details.\n\nHere is a brief description of each possible return value:\n
	CV_SUCCESS:       CVode succeeded and no roots were found.  \n\tCV_ROOT_RETURN:   CVode succeeded, and found one or more roots. If nrtfn > 1, call CVodeGetRootInfo to see which g_i were found to have a root at (*tret).  \n\tCV_TSTOP_RETURN:  CVode succeeded and returned at tstop.  \n\tCV_MEM_NULL:      The cvode_mem argument was NULL.  \n\tCV_NO_MALLOC:     cvode_mem was not allocated.  \n\tCV_ILL_INPUT:     One of the inputs to CVode is illegal. This includes the situation when a component of the error weight vectors becomes < 0 during internal time-stepping.  It also includes the situation where a root of one of the root functions was found both at t0 and very near t0.  The ILL_INPUT flag will also be returned if the linear solver routine CV--- (called by the user after calling CVodeCreate) failed to set one of the linear solver-related fields in cvode_mem or if the linear solver's init routine failed. In any case, the user should see the printed error message for more details.  \n\tCV_TOO_MUCH_WORK: The solver took mxstep internal steps but could not reach tout. The default value for mxstep is MXSTEP_DEFAULT = 500.  \n\tCV_TOO_MUCH_ACC:  The solver could not satisfy the accuracy demanded by the user for some internal step.  \n\tCV_ERR_FAILURE:   Error test failures occurred too many times (= MXNEF = 7) during one internal time step or occurred with |h| = hmin.  \n\tCV_CONV_FAILURE:  Convergence test failures occurred too many times (= MXNCF = 10) during one internal time step or occurred with |h| = hmin.  \n\tCV_LINIT_FAIL:    The linear solver's initialization function failed.  \n\tCV_LSETUP_FAIL:   The linear solver's setup routine failed in an unrecoverable manner.  \n\tCV_LSOLVE_FAIL:   The linear solver's solve routine failed in an unrecoverable manner."""
	ret = cvode.CVode(cvodememobj.obj, tout, yout.data, tret, itask)
	return ret
cvode.CVode.argtypes = [ctypes.c_void_p, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype), ctypes.c_int]
cvode.CVode.restype = ctypes.c_int

def CVodeGetDky(cvodememobj, t, k, dky):
	"""CVodeGetDky computes the kth derivative of the y function at time t, where tn-hu <= t <= tn, tn denotes the current internal time reached, and hu is the last internal step size successfully used by the solver. The user may request k=0, 1, ..., qu, where qu is the order last used. The derivative vector is returned in dky. This vector must be allocated by the caller. It is only legal to call this function after a successful return from CVode.\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate() \n\tt\t\tthe time at which the kth derivative of y is evaluated. The legal range for t is [tn-hu,tn] as described in the CVODE documentation. [float] \n\tk\t\tthe order of the derivative of y to be computed. The legal range for k is [0,qu] as described in the CVODE documentation. [int] \n\tdky\t\t the output derivative vector [((d/dy)^k)y](t). [NVector]"""
	ret = cvode.CVodeGetDky(cvodememobj.obj, t, k, dky.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetDky() failed with flag %i"%(ret))
cvode.CVodeGetDky.argtypes = [ctypes.c_void_p, realtype, ctypes.c_int, ctypes.POINTER(nvecserial._NVector)] 
cvode.CVodeGetDky.restype = ctypes.c_int

def CVodeGetWorkSpace(cvodememobj):
	"""CVodeGetWorkSpace returns the CVODE real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
	lenrw = ctypes.c_long()
	leniw = ctypes.c_long()
	ret = cvode.CVodeGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrw), ctypes.byref(leniw))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetWorkSpace() failed with flag %i"%(ret))
	return (lenrw, leniw)
cvode.CVodeGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetWorkSpace.restype = ctypes.c_int

def CVodeGetNumSteps(cvodememobj):
	"""CVodeGetNumSteps returns the cumulative number of internal steps taken by the solver\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumSteps(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumSteps() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumSteps.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumSteps.restype = ctypes.c_int

def CVodeGetNumRhsEvals(cvodememobj):
	"""CVodeGetNumRhsEvals returns the number of calls to the user's f function\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumRhsEvals.restype = ctypes.c_int

def CVodeGetNumLinSolvSetups(cvodememobj):
	"""CVodeGetNumLinSolvSetups returns the number of calls made to the linear solver's setup routine\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumLinSolvSetups(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumLinSolvSetups() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumLinSolvSetups.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumLinSolvSetups.restype = ctypes.c_int

def CVodeGetNumErrTestFails(cvodememobj):
	"""CVodeGetNumErrTestFails returns the number of local error test failures that have occured\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumErrTestFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumErrTestFails() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumErrTestFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumErrTestFails.restype = ctypes.c_int

def CVodeGetLastOrder(cvodememobj):
	"""CVodeGetLastOrder returns the order used during the last internal step\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVodeGetLastOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastOrder() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetLastOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetLastOrder.restype = ctypes.c_int

def CVodeGetCurrentOrder(cvodememobj):
	"""CVodeGetCurrentOrder returns the order to be used on the next internal step\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_int()
	ret = cvode.CVodeGetCurrentOrder(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentOrder() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentOrder.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetCurrentOrder.restype = ctypes.c_int

def CVodeGetNumStabLimOrderReds(cvodememobj):
	"""CVodeGetNumStabLimOrderReds returns the number of order reductions due to stability limit detection\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long()
	ret = cvode.CVodeGetNumStabLimOrderReds(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumStabLimOrderReds() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumStabLimOrderReds.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumStabLimOrderReds.restype = ctypes.c_int

def CVodeGetActualInitStep(cvodememobj):
	"""CVodeGetActualInitStep returns the actual initial step size used by CVODE\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetActualInitStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetActualInitStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetActualInitStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetActualInitStep.restype = ctypes.c_int

def CVodeGetLastStep(cvodememobj):
	"""CVodeGetLastStep returns the step size for the last internal step\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetLastStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetLastStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetLastStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetLastStep.restype = ctypes.c_int

def CVodeGetCurrentStep(cvodememobj):
	"""CVodeGetCurrentStep returns the step size to be attempted on the next internal step\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetCurrentStep(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentStep() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentStep.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetCurrentStep.restype = ctypes.c_int

def CVodeGetCurrentTime(cvodememobj):
	"""CVodeGetCurrentTime returns the current internal time reached by the solver\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetCurrentTime(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetCurrentTime() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetCurrentTime.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetCurrentTime.restype = ctypes.c_int

def CVodeGetTolScaleFactor(cvodememobj):
	"""CVodeGetTolScaleFactor returns a suggested factor by which the user's tolerances should be scaled when too much accuracy has been requested for some internal step\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = realtype()
	ret = cvode.CVodeGetTolScaleFactor(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetTolScaleFactor() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetTolScaleFactor.argtypes = [ctypes.c_void_p, ctypes.POINTER(realtype)]
cvode.CVodeGetTolScaleFactor.restype = ctypes.c_int

def CVodeGetErrWeights(cvodememobj, eweight):
	"""CVodeGetErrWeights returns the current error weight vector in eweight.\n\tcvodememobj\ta CvodeMemObj as returned by CVodeCreate()\n\teweight\t\tis set to the current error weight vector [NVector]"""
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
	"""CVodeGetEstLocalErrors returns the vector of estimated local errors in ele.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tele\t\tis set to the vector of current estimated local errors [NVector]"""
	ret = cvode.CVodeGetEstLocalErrors(cvodememobj.obj, ele.data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetEstLocalErrors() failed with flag %i"%(ret))
cvode.CVodeGetEstLocalErrors.argtypes = [ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector)]
cvode.CVodeGetEstLocalErrors.restype = ctypes.c_int

def CVodeGetNumGEvals(cvodememobj):
	"""CVodeGetNumGEvals returns the number of calls to the user's g function (for rootfinding)\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumGEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumGEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumGEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumGEvals.restype = ctypes.c_int

def CVodeGetRootInfo(cvodememobj, numroots):
	"""CVodeGetRootInfo returns the indices for which the function g_i was found to have a root. A list of zeroes and ones is returned, where a one in index postion i indicates a root for g_i was found.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tnumroots\t\tthe number of function for which roots should be found."""
	rootsfound = (ctypes.c_int*numroots)()
	ret = cvode.CVodeGetRootInfo(cvodememobj.obj, rootsfound)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVode() failed with flag %i"%(ret))
	return [rootsfound[i] for i in range(numroots)]
cvode.CVodeGetRootInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVodeGetRootInfo.restype = ctypes.c_int
	
def CVodeGetIntegratorStats(cvodememobj):
	"""A convenience function which returns a tuple of all available integrator stats.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\nReturn tuple is of he format:\n \n\t(nsteps, nfevals, nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur) where: \n\t\tnsteps: cumulative number of internal steps taken by the solver \n\t\tnfevals: number of rhs evaluations \n\t\tnlinsetups: number of calls to the linear setup function \n\t\tnetfails: number of local error test failures \n\t\tqlast: order used during ast internal step \n\t\tqcur: order to be used on next internal step \n\t\thinused: actual initial step size ised by CVODE \n\t\thlast: last step size used by CVODE \n\t\thcur: next step size to be used by CVODE \n\t\ttcur = current time"""
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
	"""CVodeGetNumNonlinSolvIters returns the number of nonlinear solver iterations performed.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumNonlinSolvIters(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvIters() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumNonlinSolvIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumNonlinSolvIters.restype = ctypes.c_int

def CVodeGetNumNonlinSolvConvFails(cvodememobj):
	"""CVodeGetNumNonlinSolvConvFails returns the number of nonlinear convergence failures.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	retval = ctypes.c_long(0)
	ret = cvode.CVodeGetNumNonlinSolvConvFails(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNumNonlinSolvConvFails() failed with flag %i"%(ret))
	return retval.value
cvode.CVodeGetNumNonlinSolvConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNumNonlinSolvConvFails.restype = ctypes.c_int

def CVodeGetNonlinSolvStats(cvodememobj):
	"""A convenience function that provides the nonlinear solver optional outputs in a tuple (NumNonlinSolvIters, NumNonlinsolvConvFails).\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	nniters = ctypes.c_long()
	nncfails = ctypes.c_long()
	ret = cvode.CVodeGetNonlinSolvStats(cvodememobj.obj, ctypes.byref(nniters), ctypes.byref(nncfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVodeGetNonlinSolvStats() failed with flag %i"%(ret))
	return (nniters.value, nncfails.value)
cvode.CVodeGetNonlinSolvStats.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVodeGetNonlinSolvStats.restype = ctypes.c_int

def CVodeGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVODE return flag"""
	return cvode.CVodeGetReturnFlagName(flag)
cvode.CVodeGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVodeGetReturnFlagName.restype = ctypes.c_char_p

#def CVodeFree(cvodememobj):
#	"""CVodeFree frees the problem memory cvodememobj allocated by CVodeCreate and CVodeMalloc."""
#	cvode.CVodeFree(cvodememobj.obj)
#	del cvodememobj.obj

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
	ret = cvode.ModifiedGS(v.data, h, k, p, new_vk_norm)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ModifiedGS() failed with flag %i"%(ret))
cvode.ModifiedGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype)]
cvode.ModifiedGS.restype = ctypes.c_int

def ClassicalGS(v, h, k, p, new_vk_norm, temp, s):
	ret = cvode.ClassicalGS(v.data, h, k, p, new_vk_norm, temp.data, s)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: ClassicalGS() failed with flag %i"%(ret))
cvode.ClassicalGS.argtypes = [ctypes.POINTER(ctypes.POINTER(nvecserial._NVector)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_int, ctypes.c_int, ctypes.POINTER(realtype), ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(realtype)]
cvode.ClassicalGS.restype = ctypes.c_int

def QRfact(n, h, q, job):
	ret = cvode.QRfact(n, h, q, job)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRfact() failed with flag %i"%(ret))
cvode.QRfact.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.c_int]
cvode.QRfact.restype = ctypes.c_int

def QRsol(n, h, q, b):
	ret = cvode.QRsol(n, h, q, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: QRsol() failed with flag %i"%(ret))
cvode.QRsol.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(realtype), ctypes.POINTER(realtype)]
cvode.QRsol.restype = ctypes.c_int


#########################
# sundials_smalldense.h #
#########################

def denalloc(m, n):
	ret = cvode.denalloc(m, n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denalloc could not allocate memory')
cvode.denalloc.argtypes = [ctypes.c_long, ctypes.c_long]
cvode.denalloc.restype = ctypes.POINTER(ctypes.POINTER(realtype))

def denallocpiv(n):
	ret = cvode.denallocpiv(n)
	if ret is not None:
		return ret
	else:
		raise ValueError('denallocpiv could not allocate memory')
cvode.denallocpiv.argtypes = [ctypes.c_long]
cvode.denallocpiv.restype = ctypes.POINTER(ctypes.c_long)

def denGETRF(a, m, n, p):
	return cvode.denGETRF(a, m, n, p)
cvode.denGETRF.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long, ctypes.POINTER(ctypes.c_long)]
cvode.denGETRF.restype = ctypes.c_long

def denGETRS(a, n, p, b):
	return cvode.denGETRS(a, n, p, b)
cvode.denGETRS.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype)]
cvode.denGETRS.restype = None

def denzero(a, m, n):
	cvode.denzero(a, m, n)
cvode.denzero.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.denzero.restype = None

def dencopy(a, b, m, n):
	cvode.dencopy(a, b, m, n)
cvode.dencopy.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.dencopy.restype = None

def denscale(c, a, m, n):
	cvode.denscale(c, a, m, n)
cvode.denscale.argtypes = [realtype, ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long, ctypes.c_long]
cvode.denscale.restype = None

def denaddI(a, n):
	cvode.denaddI(a, n)
cvode.denaddI.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype)), ctypes.c_long]
cvode.denaddI.restype = None

def denfreepiv(p):
	cvode.denfreepiv(p)
cvode.denfreepiv.argtypes = [ctypes.POINTER(ctypes.c_long)]
cvode.denfreepiv.restype = None

def denfree(a):
	cvode.denfree(a)
cvode.denfree.argtypes = [ctypes.POINTER(ctypes.POINTER(realtype))]
cvode.denfree.restype = None

def denprint(a, m, n):
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
	return cvode.SpbcgMalloc(l_max, vec_tmpl.data)
cvode.SpbcgMalloc.argtypes = [ctypes.c_int, ctypes.POINTER(nvecserial._NVector)]
cvode.SpbcgMalloc.restype = SpbcgMem

def SpbcgSolve(mem, A_data, x, b, pretype, delta, P_data, sx, sb, atimes, psolve, res_norm, nli, nps):
	ret = cvode.SpbcgSolve(mem, A_data, x.data, b.data, pretype, delta, P_data, sx.data, sb.data, atimes, psolve, res_norm, nli, nps)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: SpbcgSolve() failed with flag %i"%(ret))
cvode.SpbcgSolve.argtypes = [SpbcgMem, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ctypes.c_int, realtype, ctypes.c_void_p, ctypes.POINTER(nvecserial._NVector), ctypes.POINTER(nvecserial._NVector), ATimesFn, PSolveFn, ctypes.POINTER(realtype), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
cvode.SpbcgSolve.restype = ctypes.c_int

def SpbcgFree(mem):
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
	"""Links the main integrator with the CVDIAG linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	ret = cvode.CVDiag(cvodememobj.obj)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiag() failed with flag %i"%(ret))
cvode.CVDiag.argtypes = [ctypes.c_void_p]
cvode.CVDiag.restype = ctypes.c_int

def CVDiagGetWorkSpace(cvodememobj):
	"""CVDiagGetWorkSpace returns the CVDIAG real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVDiagGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVDiagGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVDiagGetWorkSpace.restype = ctypes.c_int

def CVDiagGetNumRhsEvals(cvodememobj):
	"""CVDiagGetNumRhsEvals returns the number of calls to the user f routine due to finite difference Jacobian evaluation.  Note: The number of diagonal approximate Jacobians formed is equal to the number of CVDiagSetup calls. This number is available through CVodeGetNumLinSolvSetups."""
	retval = ctypes.c_long(0)
	ret = cvode.CVDiagGetNumRhsEvals(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDiagGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVDiagGetNumRhsEvals.restype = ctypes.c_int

def CVDiagGetLastFlag(cvode_mem):
	"""CVDiagGetLastFlag returns the last error flag set by any of the CVDIAG interface functions."""
	retval = ctypes.c_int(0)
	ret = cvode.CVDiagGetLastFlag(cvode_mem, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDiagGetLastFlag() failed with flag %i"%(ret))
	return retval
cvode.CVDiagGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVDiagGetLastFlag.restype = ctypes.c_int

def CVDiagGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDIAG return flag."""
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
	"""Creates a wrapper around a python callable object, that can be used as a callback for the local approximation RHS function. Local approximation RHS functions take Nlocal (int = size of local vector), time_step (float), y (NVector), g (NVector), and f_data (c_void_p) as parameters, and return an integer."""
	if func == None:
		return ctypes.cast(None, CVLocalFn)
	exec 'def __CallbackInterface_%s(Nlocal, t, y, g, f_data):\n\treturn __ActualCallback[%i](Nlocal, t, nvecserial.NVector(y), nvecserial.NVector(g), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVLocalFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

CVCommFn = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_long, realtype, ctypes.POINTER(nvecserial._NVector), ctypes.c_void_p)
def WrapCallbackCVCommFn(func):
	"""Creates a wrapper around a python callable object, that can be used to perform all IPC necessary to approximate the RHS function. Such functions take Nlocal (int = size of local vector), time_step (float), y (NVector), and f_data (c_void_p) as parameters, and return nothing."""
	if func == None:
		return ctypes.cast(None, CVCommFn)
	exec 'def __CallbackInterface_%s(Nlocal, t, y, f_data):\n\treturn __ActualCallback[%i](Nlocal, t, nvecserial.NVector(y), f_data)'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVCommFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVBBDPrecAlloc(cvodememobj, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn):
	"""CVBBDPrecAlloc allocates and initializes a CVBBDData structure to be passed to CVSp* (and used by CVBBDPrecSetup and CVBBDPrecSolve).\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate() \n\tNlocal\t\tthe length of the local block of the vectors y etc.  on the current processor.  \n\tmudq, mldq\tthe upper and lower half-bandwidths to be used in the difference quotient computation of the local Jacobian block.  \n\tmukeep, mlkeep\tthe upper and lower half-bandwidths of the retained banded approximation to the local Jacobian block.  \n\tdqrely\t\tan optional input. It is the relative increment in components of y used in the difference quotient approximations. To specify the default, pass 0.  The default is dqrely = sqrt(unit roundoff).  \n\tgloc\t\tthe name of the user-supplied function g(t,y) that approximates f and whose local Jacobian blocks are to form the preconditioner.  \n\tcfn\t\tthe name of the user-defined function that performs necessary interprocess communication for the execution of gloc.  \n\tCVBBDPrecAlloc returns an object reperesenting the storage allocated or None if the request for storage cannot be satisfied."""
	ret = cvode.CVBBDPrecAlloc(cvodememobj.obj, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, WrapCallbackCVLocalFn(gloc), WrapCallbackCVCommFn(cfn))
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBBDPrecAlloc() failed to allocate memory")
	return ret
cvode.CVBBDPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, ctypes.c_long, realtype, CVLocalFn, CVCommFn]
cvode.CVBBDPrecAlloc.restype = ctypes.c_void_p

def CVBBDSptfqmr(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSptfqmr links the CVBBDPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions: \n\t1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory \n\t2) Sets the preconditioner data structure for CVSPTFQMR \n\t3) Sets the preconditioner setup routine for CVSPTFQMR \n\t4) Sets the preconditioner solve routine for CVSPTFQMR \n\nIts first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSptfqmr.  \n\nPossible return values are: \n\tCVSPILS_SUCCESS      if successful \n\tCVSPILS_MEM_NULL     if the cvode memory was NULL \n\tCVSPILS_LMEM_NULL    if the cvsptfqmr memory was NULL \n\tCVSPILS_MEM_FAIL     if there was a memory allocation failure \n\tCVSPILS_ILL_INPUT    if a required vector operation is missing \n\tCVBBDPRE_PDATA_NULL  if the bbd_data was NULL"""
	ret = cvode.CVBBDSptfqmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSptfqmr() failed with flag %i"%(ret))
cvode.CVBBDSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSptfqmr.restype = ctypes.c_int

def CVBBDSpbcg(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSptfqmr links the CVBBDPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions:\n\t1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory;\n\t2) Sets the preconditioner data structure for CVSPTFQMR\n\t3) Sets the preconditioner setup routine for CVSPTFQMR\n\t4) Sets the preconditioner solve routine for CVSPTFQMR\n\nIts first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSptfqmr."""
	ret = cvode.CVBBDSpbcg(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpbcg() failed with flag %i"%(ret))
cvode.CVBBDSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSpbcg.restype = ctypes.c_int

def CVBBDSpgmr(cvodememobj, pretype, maxl, bbd_data):
	"""CVBBDSpbcg links the CVBBDPRE preconditioner to the CVSPBCG linear solver. It performs the following actions:\n\t1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory;\n\t2) Sets the preconditioner data structure for CVSPBCG\n\t3) Sets the preconditioner setup routine for CVSPBCG\n\t4) Sets the preconditioner solve routine for CVSPBCG\n\nIts first 3 arguments are the same as for CVSpbcg (see cvspbcg.h). The last argument is the pointer to the CVBBDPRE memory block returned by CVBBDPrecAlloc. Note that the user need not call CVSpbcg."""
	ret = cvode.CVBBDSpgmr(cvodememobj.obj, pretype, maxl, bbd_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBBDSpgmr() failed with flag %i"%(ret))
cvode.CVBBDSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBBDSpgmr.restype = ctypes.c_int

def CVBBDPrecReInit(bbd_data, mudq, mldq, dqrely, gloc, cfn):
	"""CVBBDPrecReInit re-initializes the BBDPRE module when solving a sequence of problems of the same size with CVSPGMR/CVBBDPRE or CVSPBCG/CVBBDPRE or CVSPTFQMR/CVBBDPRE provided there is no change in Nlocal, mukeep, or mlkeep. After solving one problem, and after calling CVodeReInit to re-initialize the integrator for a subsequent problem, call CVBBDPrecReInit. Then call CVSpgmrSet* or CVSpbcgSet* or CVSptfqmrSet* functions if necessary for any changes to CVSPGMR, CVSPBCG, or CVSPTFQMR parameters, before calling CVode.\n\nThe first argument to CVBBDPrecReInit must be the pointer pdata that was returned by CVBBDPrecAlloc. All other arguments have the same names and meanings as those of CVBBDPrecAlloc."""
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
	"""CVBBDPrecGetWorkSpace returns the CVBBDPrec real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
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
	"""Returns the name of the constant associated with a CVBBDPRE return flag."""
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
	"""CVBandPrecAlloc allocates and initializes a CVBandPrecData structure to be passed to CVSp* (and subsequently used by CVBandPrecSetup and CVBandPrecSolve). The parameters of CVBandPrecAlloc are as follows:\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate() \n\tN\tthe problem size.  \n\tmu\tthe upper half bandwidth.  \n\tml\tthe lower half bandwidth.\nCVBandPrecAlloc returns the storage pointer of type CVBandPrecData, or NULL if the request for storage cannot be satisfied.  \nNOTE:\tThe band preconditioner assumes a serial implementation of the NVECTOR package. Therefore, CVBandPrecAlloc will first test for a compatible N_Vector internal representation by checking for required functions."""
	ret = cvode.CVBandPrecAlloc(cvodememobj.obj, N, mu, ml)
	if ret is None:
		raise AssertionError("SUNDIALS ERROR: CVBandPrecAlloc() failed to allocate memory")
	return ret
cvode.CVBandPrecAlloc.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvode.CVBandPrecAlloc.restype = ctypes.c_void_p

def CVBPSptfqmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSptfqmr links the CVBANDPPRE preconditioner to the CVSPTFQMR linear solver. It performs the following actions: \n\t1) Calls the CVSPTFQMR specification routine and attaches the CVSPTFQMR linear solver to the integrator memory; \n\t2) Sets the preconditioner data structure for CVSPTFQMR \n\t3) Sets the preconditioner setup routine for CVSPTFQMR \n\t4) Sets the preconditioner solve routine for CVSPTFQMR \nIts first 3 arguments are the same as for CVSptfqmr (see cvsptfqmr.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSptfqmr."""
	ret = cvode.CVBPSptfqmr(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSptfqmr() failed with flag %i"%(ret))
cvode.CVBPSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBPSptfqmr.restype = ctypes.c_int

def CVBPSpbcg(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpbcg links the CVBANDPPRE preconditioner to the CVSPBCG linear solver. It performs the following actions: \n\t1) Calls the CVSPBCG specification routine and attaches the CVSPBCG linear solver to the integrator memory; \n\t2) Sets the preconditioner data structure for CVSPBCG \n\t3) Sets the preconditioner setup routine for CVSPBCG \n\t4) Sets the preconditioner solve routine for CVSPBCG \nIts first 3 arguments are the same as for CVSpbcg (see cvspbcg.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSpbcg."""
	ret = cvode.CVBPSpbcg(cvodememobj.obj, pretype, maxl, p_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBPSpbcg() failed with flag %i"%(ret))
cvode.CVBPSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
cvode.CVBPSpbcg.restype = ctypes.c_int

def CVBPSpgmr(cvodememobj, pretype, maxl, p_data):
	"""CVBPSpgmr links the CVBANDPPRE preconditioner to the CVSPGMR linear solver. It performs the following actions: \n\t1) Calls the CVSPGMR specification routine and attaches the CVSPGMR linear solver to the integrator memory; \n\t2) Sets the preconditioner data structure for CVSPGMR \n\t3) Sets the preconditioner setup routine for CVSPGMR \n\t4) Sets the preconditioner solve routine for CVSPGMR \n\tIts first 3 arguments are the same as for CVSpgmr (see cvspgmr.h). The last argument is the pointer to the CVBANDPPRE memory block returned by CVBandPrecAlloc. Note that the user need not call CVSpgmr."""
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
	"""CVBandPrecGetWorkSpace returns the CVBandPrec real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
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
	"""Returns the name of the constant associated with a CVBBDPRE return flag."""
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
		"""Instantiates a new dense representation of a banded matrix.\n\tinit\tdimension of matrix to create; always square [int]\n\tmu\twidth of portion of band above diagonal\n\tml\twidth of portion of band below diagonal\n\tsmu\tis the storage upper bandwidth, mu <= smu <= size-1.  The BandGBTRF routine writes the LU factors into the storage for A. The upper triangular factor U, however, may have an upper bandwidth as big as MIN(size-1,mu+ml) because of partial pivoting. The smu field holds the upper bandwidth allocated for A."""
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
	"""BandAllocPiv allocates memory for pivot information to be filled in by the BandGBTRF routine during the factorization of an N by N band matrix. Returns a pointer, which should be passed as is to the [CV]BandPiv* functions.\n\tN\tthe size of the banded matrix for which to allocate pivot information space [int]"""
	return cvode.BandAllocPiv(N)
cvode.BandAllocPiv.argtypes = [ctypes.c_int]
cvode.BandAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def BandGBTRF(A, p):
	"""BandGBTRF performs the LU factorization of the N by N band matrix A. This is done using standard Gaussian elimination with partial pivoting.\n\nA successful LU factorization leaves the "matrix" A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower triangular matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower triangular part of A contains the multipliers, I-L.\n\nBandGBTRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero.\n\nImportant Note: A must be allocated to accommodate the increase in upper bandwidth that occurs during factorization. If mathematically, A is a band matrix with upper bandwidth mu and lower bandwidth ml, then the upper triangular factor U can have upper bandwidth as big as smu = MIN(n-1,mu+ml). The lower triangular factor L has lower bandwidth ml. Allocate A with call A = BandAllocMat(N,mu,ml,smu), where mu, ml, and smu are as defined above. The user does not have to zero the "extra" storage allocated for the purpose of factorization. This will handled by the BandGBTRF routine.  ret = cvode.BandGBTRF(A, p)"""
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
	"""Zeroes the entire matrix\n\tA\tthe matrix [BandMat]"""
	cvode.BandZero(A.data)
cvode.BandZero.argtypes = [_BandMat] 
cvode.BandZero.restype = None

def BandCopy(A, B, copymu, copyml):
	"""Copies the contents of banded matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [BandMat]\n\tB\tthe source matrix [BandMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	cvode.BandCopy(A.data, B.data, copymu, copyml)
cvode.BandCopy.argtypes = [_BandMat, _BandMat, ctypes.c_long, ctypes.c_long]
cvode.BandCopy.restype = None

def BandScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [BandMat]\n\tc\tthe scale factor [float]"""
	cvode.BandScale(c, A.data)
cvode.BandScale.argtypes = [realtype, _BandMat]
cvode.BandScale.restype = None

def BandAddI(A):
	"""Adds 1.0 to the diagonal of the banded matrix A\n\tA\tthe matrix [BandMat]"""
	cvode.BandAddI(A.data)
cvode.BandAddI.argtypes = [_BandMat]
cvode.BandAddI.restype = None

def BandFreeMat(A):
	"""BandFreeMat frees the memory allocated by BandAllocMat for the band matrix A.\n\tA\tthe matrix [BandMat]"""
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
	"""A call to the CVBand function links the main CVODE integrator with the CVBAND linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate() \n\tN\t\tthe size of the ODE system.  \n\tmupper\t\tthe upper bandwidth of the band Jacobian approximation.\n\tmlower\t\tthe lower bandwidth of the band Jacobian approximation."""
	ret = cvode.CVBand(cvodememobj.obj, N, mupper, mlower)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvode.CVBand.argtypes = [ctypes.c_void_p, ctypes.c_long, ctypes.c_long, ctypes.c_long]
cvode.CVBand.restype = ctypes.c_int

def CVBandSetJacFn(cvodememobj, func, jac_data):
	"""CVBandSetJacFn sets the band Jacobian approximation function to be used.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tfunc\t\ta python callable suitable for wrapping by WrapCallbackCVBandJacFn\n\tjac_data\ta pointer to user data"""
	ret = cvode.CVBandSetJacFn(cvodememobj.obj, WrapCallbackCVBandJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBand() failed with flag %i"%(ret))
cvode.CVBandSetJacFn.argtypes = [ctypes.c_void_p, CVBandJacFn, ctypes.c_void_p]
cvode.CVBand.restype = ctypes.c_int

def CVBandGetWorkSpace(cvodememobj):
	"""CVBandGetWorkSpace returns the CVBand real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
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
	"""CVBandGetLastFlag returns the last error flag set by any of the CVBAND interface functions."""
	retval = ctypes.c_int()
	ret = cvode.CVBandGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVBandGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVBandGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVBandGetLastFlag.restype = ctypes.c_int

def CVBandGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVBAND return flag."""
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
		"""Instantiates a dense matrix object\n\tM\tnumber of rows\n\tN\tnumber of columns"""
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
	"""DenseAllocPiv allocates memory for pivot information to be filled in by the DenseGETRF routine during the factorization of an N by N dense matrix. Returns a pointer, which should be passed as is to the [CV]DensePiv* functions.\n\tN\tthe size of the dense matrix for which to allocate pivot information space [int]"""
	return cvode.DenseAllocPiv(N)
cvode.DenseAllocPiv.argtypes = [ctypes.c_int]
cvode.DenseAllocPiv.restype = ctypes.POINTER(ctypes.c_long)

def DenseGETRF(A, p):
	"""DenseGETRF performs the LU factorization of the M by N dense matrix A. This is done using standard Gaussian elimination with partial (row) pivoting. Note that this only applies to matrices with M >= N and full column rank.\n\nA successful LU factorization leaves the matrix A and the pivot array p with the following information:\n\t1 p[k] contains the row number of the pivot element chosen at the beginning of elimination step k, k=0, 1, ..., N-1.\n\t2 If the unique LU factorization of A is given by PA = LU, where P is a permutation matrix, L is a lower trapezoidal matrix with all 1's on the diagonal, and U is an upper triangular matrix, then the upper triangular part of A (including its diagonal) contains U and the strictly lower trapezoidal part of A contains the multipliers, I-L.\n\nFor square matrices (M=N), L is unit lower triangular.\n\nDenseGETRF returns 0 if successful. Otherwise it encountered a zero diagonal element during the factorization. In this case it returns the column index (numbered from one) at which it encountered the zero."""
	ret = cvode.DenseGETRF(A.data, p)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRF() failed with flag %i"%(ret))
cvode.DenseGETRF.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long)]
cvode.DenseGETRF.restype = ctypes.c_long

def DenseGETRS(A, p, b):
	"""DenseGETRS solves the N-dimensional system A x = b using the LU factorization in A and the pivot information in p computed in DenseGETRF. The solution x is returned in b. This routine cannot fail if the corresponding call to DenseGETRF did not fail.\nDenseGETRS does NOT check for a squre matrix!"""
	ret = cvode.DenseGETRS(A.data, p, b)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: DenseGETRS() failed with flag %i"%(ret))
cvode.DenseGETRS.argtypes = [_DenseMat, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(realtype )]
cvode.DenseGETRS.restype = None

def DenseZero(A):
	"""Zeroes the entire matrix\n\tA\tthe matrix [DenseMat]"""
	cvode.DenseZero(A.data)
cvode.DenseZero.argtypes = [_DenseMat] 
cvode.DenseZero.restype = None

def DenseCopy(A, B):
	"""Copies the contents of dense matrix B into A, overwriting A's contents.\n\tA\tthe destination matrix [DenseMat]\n\tB\tthe source matrix [DenseMat]\n\tcopymu\tupper band width of matrices [int]\n\tcopyml\tlower band width of matrices [int]"""
	cvode.DenseCopy(A.data, B.data)
cvode.DenseCopy.argtypes = [_DenseMat, _DenseMat]
cvode.DenseCopy.restype = None

def DenseScale(c, A):
	"""Scales the matrix A by c\n\tA\tthe matrix [DenseMat]\n\tc\tthe scale factor [float]"""
	cvode.DenseScale(c, A.data)
cvode.DenseScale.argtypes = [realtype, _DenseMat]
cvode.DenseScale.restype = None

def DenseAddI(A):
	"""Adds 1.0 to the diagonal of the dense matrix A\n\tA\tthe matrix [DenseMat]"""
	cvode.DenseAddI(A.data)
cvode.DenseAddI.argtypes = [_DenseMat]
cvode.DenseAddI.restype = None

def DenseFreeMat(A):
	"""DenseFreeMat frees the memory allocated by DenseAllocMat for the band matrix A.\n\tA\tthe matrix [DenseMat]"""
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
	if (func == None):
		return ctypes.cast(None, CVDenseJacFn)
	exec 'def __CallbackInterface_%s(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):\n\treturn __ActualCallback[%i](N, DenseMat(J), t, nvecserial.NVector(y), nvecserial.NVector(fy), jac_data, nvecserial.NVector(tmp1), nvecserial.NVector(tmp2), nvecserial.NVector(tmp3))'%(func.func_name, len(__ActualCallback))
	__ActualCallback.append(func)
	tmp = CVDenseJacFn(eval("__CallbackInterface_%s"%(func.func_name)))
	__Callback.append(tmp)
	return tmp

def CVDense(cvodememobj, N):
	"""A call to the CVDense function links the main CVODE integrator with the CVDENSE linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tN\t\tthe size of the ODE system."""
	ret = cvode.CVDense(cvodememobj.obj, N)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvode.CVDense.argtypes = [ctypes.c_void_p, ctypes.c_long]
cvode.CVDense.restype = ctypes.c_int

def CVDenseSetJacFn(cvodememobj, func, jac_data):
	"""CVDenseSetJacFn sets the dense Jacobian approximation function to be used.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tfunc\t\ta python callable suitable for wrapping by WrapCallbackCVDenseJacFn\n\tjac_data\ta pointer to user data"""
	ret = cvode.CVDenseSetJacFn(cvodememobj.obj, WrapCallbackCVDenseJacFn(func), jac_data)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDense() failed with flag %i"%(ret))
cvode.CVDenseSetJacFn.argtypes = [ctypes.c_void_p, CVDenseJacFn, ctypes.c_void_p]
cvode.CVDense.restype = ctypes.c_int

def CVDenseGetWorkSpace(cvodememobj):
	"""CVDenseGetWorkSpace returns the CVDense real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
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
	"""CVDenseGetLastFlag returns the last error flag set by any of the CVDENSE interface functions."""
	retval = ctypes.c_int()
	ret = cvode.CVDenseGetLastFlag(cvodememobj.obj, ctypes.byref(retval))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVDenseGetNumRhsEvals() failed with flag %i"%(ret))
	return retval.value
cvode.CVDenseGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVDenseGetLastFlag.restype = ctypes.c_int

def CVDenseGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVDENSE return flag."""
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
	"""CVSpilsGetWorkSpace returns the CVSpils real and integer workspaces as a tuple (integer, real)\n\tcvodememobj\ta CVodeMemObj as returned by CVodeCreate()"""
	lenrwLS = ctypes.c_long(0)
	leniwLS = ctypes.c_long(0)
	ret = cvode.CVSpilsGetWorkSpace(cvodememobj.obj, ctypes.byref(lenrwLS), ctypes.byref(leniwLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetWorkSpace() failed with flag %i"%(ret))
	return (lenrwLS.value, leniwLS.value)
cvode.CVSpilsGetWorkSpace.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetWorkSpace.restype = ctypes.c_int

def CVSpilsGetNumPrecEvals(cvodememobj):
	"""CVSpilsGetNumPrecEvals returns the number of preconditioner evaluations, i.e. the number of calls made to PrecSetup with jok==FALSE.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	npevals = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumPrecEvals(cvodememobj.obj, ctypes.byref(npevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecEvals() failed with flag %i"%(ret))
	return npevals.value
cvode.CVSpilsGetNumPrecEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumPrecEvals.restype = ctypes.c_int

def CVSpilsGetNumPrecSolves(cvodememobj):
	"""CVSpilsGetNumPrecSolves returns the number of calls made to PrecSolve.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	npsolves = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumPrecSolves(cvodememobj.obj, ctypes.byref(npsolves))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumPrecSolves() failed with flag %i"%(ret))
	return npsolves.value
cvode.CVSpilsGetNumPrecSolves.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumPrecSolves.restype = ctypes.c_int

def CVSpilsGetNumLinIters(cvodememobj):
	"""CVSpilsGetNumLinIters returns the number of linear iterations.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	nliters = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumLinIters(cvodememobj.obj, ctypes.byref(nliters))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumLinIters() failed with flag %i"%(ret))
	return nliters.value
cvode.CVSpilsGetNumLinIters.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumLinIters.restype = ctypes.c_int

def CVSpilsGetNumConvFails(cvodememobj):
	"""CVSpilsGetNumConvFails returns the number of linear convergence failures.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	nlcfails = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumConvFails(cvodememobj.obj, ctypes.byref(nlcfails))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumConvFails() failed with flag %i"%(ret))
	return nlcfails.value
cvode.CVSpilsGetNumConvFails.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumConvFails.restype = ctypes.c_int

def CVSpilsGetNumJtimesEvals(cvodememobj):
	"""CVSpilsGetNumJtimesEvals returns the number of calls to jtimes.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	njvevals = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumJtimesEvals(cvodememobj.obj, ctypes.byref(njvevals))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumJtimesEvals() failed with flag %i"%(ret))
	return njvevals.value
cvode.CVSpilsGetNumJtimesEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumJtimesEvals.restype = ctypes.c_int

def CVSpilsGetNumRhsEvals(cvodememobj):
	"""CVSpilsGetNumRhsEvals returns the number of calls to the user f routine due to finite difference Jacobian times vector evaluation.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()"""
	nfevalsLS = ctypes.c_long(0)
	ret = cvode.CVSpilsGetNumRhsEvals(cvodememobj.obj, ctypes.byref(nfevalsLS))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetNumRhsEvals() failed with flag %i"%(ret))
	return nfevalsLS.value
cvode.CVSpilsGetNumRhsEvals.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_long)]
cvode.CVSpilsGetNumRhsEvals.restype = ctypes.c_int

def CVSpilsGetLastFlag(cvodememobj):
	"""CVSpilsGetLastFlag returns the last error flag set by any of the CVSPILS interface functions."""
	flag = ctypes.c_int(0)
	ret = cvode.CVSpilsGetLastFlag(cvodememobj.obj, ctypes.byref(flag))
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpilsGetLastFlag() failed with flag %i"%(ret))
	return flag.value
cvode.CVSpilsGetLastFlag.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
cvode.CVSpilsGetLastFlag.restype = ctypes.c_int

def CVSpilsGetReturnFlagName(flag):
	"""Returns the name of the constant associated with a CVSPILS return flag."""
	return cvode.CVSpilsGetReturnFlagName(flag)
cvode.CVSpilsGetReturnFlagName.argtypes = [ctypes.c_int]
cvode.CVSpilsGetReturnFlagName.restype = ctypes.c_char_p

##################
# cvode_spbcgs.h #
##################

def CVSpbcg(cvodememobj, pretype, maxl):
	"""A call to the CVSpbcg function links the main CVODE integrator with the CVSPBCG linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tpretype\t\tis the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in iterative.h. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.\n\tmaxl\t\tis the maximum Krylov dimension. This is an optional input to the CVSPBCG solver. Pass 0 to use the default value CVSPBCG_MAXL=5."""
	ret = cvode.CVSpbcg(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpbcg() failed with flag %i"%(ret))
cvode.CVSpbcg.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSpbcg.restype = ctypes.c_int

#################
# cvode_spgmr.h #
#################

def CVSpgmr(cvodememobj, pretype, maxl):
	"""A call to the CVSpgmr function links the main CVODE integrator with the CVSPGMR linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tpretype\t\tis the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in sundials_iterative.h.  These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.\n\tmaxl\t\tis the maximum Krylov dimension. This is an optional input to the CVSPGMR solver. Pass 0 to use the default value CVSPGMR_MAXL=5."""
	ret = cvode.CVSpgmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSpgmr() failed with flag %i"%(ret))
cvode.CVSpgmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSpgmr.restype = ctypes.c_int

###################
# cvode_sptfqmr.h #
###################

def CVSptfqmr(cvodememobj, pretype, maxl):
	"""A call to the CVSptfqmr function links the main CVODE integrator with the CVSPTFQMR linear solver.\n\tcvodememobj\ta CvodeMemObj as returned by CvodeCreate()\n\tpretype\t\tis the type of user preconditioning to be done.  This must be one of the four enumeration constants PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in iterative.h. These correspond to no preconditioning, left preconditioning only, right preconditioning only, and both left and right preconditioning, respectively.\n\tmaxl\t\tis the maximum Krylov dimension. This is an optional input to the CVSPTFQMR solver. Pass 0 to use the default value CVSPILS_MAXL=5."""
	ret = cvode.CVSptfqmr(cvodememobj.obj, pretype, maxl)
	if ret < 0:
		raise AssertionError("SUNDIALS ERROR: CVSptfqmr() failed with flag %i"%(ret))
cvode.CVSptfqmr.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
cvode.CVSptfqmr.restype = ctypes.c_int
