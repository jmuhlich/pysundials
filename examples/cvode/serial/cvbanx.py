try:
	from pysundials import cvode
except ImportError:
	import cvode
import ctypes
import math

class userdata(ctypes.Structure):
	_fields_ = [("dx", cvode.realtype), ("dy", cvode.realtype), ("hdcoef", cvode.realtype), ("hacoef", cvode.realtype), ("vdcoef", cvode.realtype)]

def f(t, u, udot, f_data):
	#print "u:"
	#for s in u:
	#	print "%11.8lf"%s
	#print
	data = ctypes.cast(f_data, ctypes.POINTER(userdata))
	hordc = data.contents.hdcoef
	horac = data.contents.hacoef
	verdc = data.contents.vdcoef
	#print hordc, horac, verdc
	
	for j in range(5):
		for i in range(10):
			uij = u[j+(i*5)]

			if j == 0:
				udn = 0.0
			else:
				udn = u[(j-1)+(i*5)]

			if j == 4:
				uup = 0.0
			else:
				uup = u[(j+1)+(i*5)]

			if i == 0:
				ult = 0.0
			else:
				ult = u[j+((i-1)*5)]

			if i == 9:
				urt = 0.0
			else:
				urt = u[j+((i+1)*5)]
			
			hdiff = hordc*(ult - 2.0*uij + urt)
			hadv = horac*(urt - ult)
			vdiff = verdc*(uup - 2.0*uij + udn)
			udot[j+(i*5)] = hdiff + hadv + vdiff
	
	#print "udot:"
	#for s in udot:
	#	print "%11.8lf"%s
	#print
	return 0

def Jac(N, mu, ml, J, t, u, fu, jac_data, tmp1, tmp2, tmp3):
	data = ctypes.cast(jac_data, ctypes.POINTER(userdata))
	hordc = data.contents.hdcoef
	horac = data.contents.hacoef
	verdc = data.contents.vdcoef

	for j in range(5):
		for i in range(10):
			k = j + i*5
			J[k][k] = -2.0*(verdc+hordc) 
			if i != 0:
				J[k][k-5]= hordc + horac
			
			if i != 9:
				J[k][k+5] = hordc - horac
			
			if j != 0:
				J[k][k-1] = verdc
			
			if j != 4:
				J[k][k+1] = verdc

	return 0

#begin main code
u = [0]*50

for j in range(5):
	y = (j+1)*1.0/(5+1)
	for i in range(10):
		x = (i+1)*2.0/(10+1)
		u[j+(i*5)] = x*(2.0 - x) * y*(1.0 - y)*math.exp(5.0*x*y)

u = cvode.NVector(u)

reltol = 0
abstol = 1.0e-5

data = userdata()
data.dx = 2.0/(10+1)
data.dy = 1.0/(5+1)
data.hdcoef = 1.0/(data.dx**2)
data.hacoef = 0.5/(2.0*data.dx)
data.vdcoef = 1.0/(data.dy**2)

cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON);
#for s in u:
#	print "%11.8lf"%s
cvode.CVodeMalloc(cvode_mem, f, 0.0, u, cvode.CV_SS, reltol, abstol)
cvode.CVodeSetFdata(cvode_mem, ctypes.pointer(data));
cvode.CVBand(cvode_mem, 50, 5, 5);
cvode.CVBandSetJacFn(cvode_mem, Jac, ctypes.pointer(data));

umax = max(abs(u))
print "\n2-D Advection-Diffusion Equation"
print "Mesh dimensions = %d X %d"%(10,5)
print "Total system size = %d"%(50)
print "Tolerance parameters: reltol = %lg   abstol = %lg\n"%(reltol, abstol)
print "At t = %lg      max.norm(u) =%14.6le"%(0.0, umax)

t = cvode.realtype(0)
iout = 1
tout = 0.1
while iout <= 10:
	cvode.CVode(cvode_mem, tout, u, ctypes.byref(t), cvode.CV_NORMAL)
	umax = max(abs(u))
	nst = cvode.CVodeGetNumSteps(cvode_mem)
  	print "At t = %4.2f   max.norm(u) =%14.6le   nst = %4ld"%(t.value, umax, nst)
	iout +=1
	tout += 0.1

nst = cvode.CVodeGetNumSteps(cvode_mem)
nfe = cvode.CVodeGetNumRhsEvals(cvode_mem)
nsetups = cvode.CVodeGetNumLinSolvSetups(cvode_mem)
netf = cvode.CVodeGetNumErrTestFails(cvode_mem)
nni = cvode.CVodeGetNumNonlinSolvIters(cvode_mem)
ncfn = cvode.CVodeGetNumNonlinSolvConvFails(cvode_mem)
nfeLS = cvode.CVBandGetNumRhsEvals(cvode_mem)
nje = cvode.CVBandGetNumJacEvals(cvode_mem)

print "\nFinal Statistics:"
print "nst = %-6i nfe  = %-6i nsetups = %-6i nfeLS = %-6i nje = %i"%(nst, nfe, nsetups, nfeLS, nje)
print "nni = %-6ld ncfn = %-6ld netf = %-6ld\n"%(nni, ncfn, netf)
