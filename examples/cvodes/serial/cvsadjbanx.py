try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math

XMAX = 2.0
YMAX = 1.0
MX = 40
MY = 20
NEQ = MX*MY
ATOL = 1.e-5
RTOLB = 1.e-6
T0 = 0.0
T1 = 0.1
DTOUT = 0.1
NOUT = 10
TOUT = 1.0
NSTEP = 50

class UserData(ctypes.Structure):
	_fields_ = [
		("dx", cvodes.realtype),
		("dy", cvodes.realtype),
		("hdcoef", cvodes.realtype),
		("hacoef", cvodes.realtype),
		("vdcoef", cvodes.realtype)
	]
PUserData = ctypes.POINTER(UserData)

def f(t, u, udot, f_data):
	data = ctypes.cast(f_data, PUserData).contents
	hordc = data.hdcoef
	horac = data.hacoef
	verdc = data.vdcoef
	
	for j in range(MY):
		for i in range(MX):
			uij = u[j+(i*MY)]

			if j == 0:
				udn = 0.0
			else:
				udn = u[(j-1)+(i*MY)]

			if j == MY-1:
				uup = 0.0
			else:
				uup = u[(j+1)+(i*MY)]

			if i == 0:
				ult = 0.0
			else:
				ult = u[j+((i-1)*MY)]

			if i == MX-1:
				urt = 0.0
			else:
				urt = u[j+((i+1)*MY)]
			
			hdiff = hordc*(ult - 2.0*uij + urt)
			hadv = horac*(urt - ult)
			vdiff = verdc*(uup - 2.0*uij + udn)
			udot[j+(i*MY)] = hdiff + hadv + vdiff
	
	return 0

def Jac(N, mu, ml, J, t, u, fu, jac_data, tmp1, tmp2, tmp3):
	data = ctypes.cast(jac_data, PUserData).contents
	hordc = data.hdcoef
	horac = data.hacoef
	verdc = data.vdcoef

	for j in range(MY):
		for i in range(MX):
			k = j + i*MY
			J[k][k] = -2.0*(verdc+hordc) 
			if i != 0:
				J[k-MY][k]= hordc + horac
			
			if i != MX-1:
				J[k+MY][k] = hordc - horac
			
			if j != 0:
				J[k-1][k] = verdc
			
			if j != MY-1:
				J[k+1][k] = verdc

	return 0

def fB(tB, u, uB, uBdot, f_dataB):
	data = ctypes.cast(f_dataB, PUserData).contents
	hordc = data.hdcoef
	horac = data.hacoef
	verdc = data.vdcoef
	
	for j in range(MY):
		for i in range(MX):
			uBij = uB[j+(i*MY)]

			if j == 0:
				uBdn = 0.0
			else:
				uBdn = uB[(j-1)+(i*MY)]

			if j == MY-1:
				uBup = 0.0
			else:
				uBup = uB[(j+1)+(i*MY)]

			if i == 0:
				uBlt = 0.0
			else:
				uBlt = uB[j+((i-1)*MY)]

			if i == MX-1:
				uBrt = 0.0
			else:
				uBrt = uB[j+((i+1)*MY)]
			
			hdiffB = hordc*(- uBlt + 2.0*uBij - uBrt)
			hadvB = horac*(uBrt - uBlt)
			vdiffB = verdc*(- uBup + 2.0*uBij - uBdn)
			uBdot[j+(i*MY)] = hdiffB + hadvB + vdiffB - 1.0
	
	return 0

def JacB(NB, muB, mlB, JB, tB, u, uB, fuB, jac_dataB, tmp1B, tmp2B, tmp3B):
	data = ctypes.cast(jac_dataB, PUserData).contents
	hordc = data.hdcoef
	horac = data.hacoef
	verdc = data.vdcoef

	for j in range(MY):
		for i in range(MX):
			k = j + i*MY
			JB[k][k] = 2.0*(verdc+hordc) 
			if i != 0:
				JB[k-MY][k]= - hordc + horac
			
			if i != MX-1:
				JB[k+MY][k] = - hordc - horac
			
			if j != 0:
				JB[k-1][k] = - verdc
			
			if j != MY-1:
				JB[k+1][k] = - verdc

	return 0

def PrintOutput(uB, data):
	x = y = 0

	dx = data.dx
	dy = data.dy

	uBmax = 0
	for j in range(MY):
		for i in range(MX):
			uBij = uB[j+(i*MY)]
			if abs(uBij) > uBmax:
				uBmax = uBij
				x = i*dx
				y = j*dy

	print "\nMaximum sensitivity"
	print "  lambda max = %le"%(uBmax)
	print "at"
	print "  x = %le\n  y = %le"%(x, y)

#begin main code
t = cvodes.realtype(0)
ncheck = ctypes.c_int(0)

reltol = 0
abstol = 1.0e-5

data = UserData()
data.dx = XMAX/(MX+1)
data.dy = YMAX/(MY+1)
data.hdcoef = 1.0/(data.dx**2)
data.hacoef = 1.5/(2.0*data.dx)
data.vdcoef = 1.0/(data.dy**2)

u = [0]*NEQ

for j in range(MY):
	y = (j+1)*data.dy
	for i in range(MX):
		x = (i+1)*data.dx
		u[j+(i*MY)] = x*(XMAX - x) * y*(YMAX - y)*math.exp(5.0*x*y)

u = cvodes.NVector(u)

print "\nCreate and allocate CVODES memory for forward runs"

cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVodeMalloc(cvode_mem, f, T0, u, cvodes.CV_SS, reltol, abstol)
cvodes.CVBand(cvode_mem, NEQ, MY, MY)
cvodes.CVBandSetJacFn(cvode_mem, Jac, ctypes.pointer(data))

print "\nAllocate global memory"

cvadj_mem = cvodes.CVadjMalloc(cvode_mem, NSTEP, cvodes.CV_HERMITE)

print "\nForward integration"
cvodes.CVodeF(cvadj_mem, TOUT, u, ctypes.byref(t), cvodes.CV_NORMAL, ctypes.byref(ncheck))

reltolB = RTOLB
abstolB = ATOL

uB = cvodes.NVector([0]*NEQ)

print "\nCreate and allocate CVODES memory for backward run"

cvodes.CVodeCreateB(cvadj_mem, cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeSetFdataB(cvadj_mem, ctypes.pointer(data))
cvodes.CVodeMallocB(cvadj_mem, fB, TOUT, uB, cvodes.CV_SS, reltolB, abstolB)
cvodes.CVBandB(cvadj_mem, NEQ, MY, MY)
cvodes.CVBandSetJacFnB(cvadj_mem, JacB, ctypes.pointer(data))

print "\nBackward integration"
cvodes.CVodeB(cvadj_mem, T0, uB, ctypes.byref(t), cvodes.CV_NORMAL)

PrintOutput(uB, data)
