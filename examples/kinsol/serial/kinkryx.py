try:
	from pysundials import kinsol
	from pysundials import nvecserial
except ImportError:
	import kinsol
	import nvecserial
import ctypes
import math

NUM_SPECIES = 6
MX = 8
MY = 8
NSMX = NUM_SPECIES * MX
NEQ = NSMX * MY
AA = 1.0
EE = 10000.0
GG = 0.5e-6
BB = 1.0
DPREY = 1.0
DPRED = 0.5
ALPHA = 1.0
AX = 1.0
AY = 1.0
FTOL = 1.e-7
STOL = 1.e-13
THOUSAND = 1000.0
PREYIN = 1.0
PREDIN = 30000.0

#define IJ_Vptr(vv,i,j)   (&NV_Ith_S(vv, i*NUM_SPECIES + j*NSMX))

class UserData(ctypes.Structure):
	_fields_ = [
		('P', (ctypes.POINTER(ctypes.POINTER(kinsol.realtype))*MX)*MY),
		('pivot', (ctypes.POINTER(ctypes.c_long)*MX)*MY),
		('acoef', ctypes.POINTER(ctypes.POINTER(kinsol.realtype))),
		('bcoef', ctypes.POINTER(kinsol.realtype)),
		('rates', nvecserial.PVector),
		('cox', ctypes.POINTER(kinsol.realtype)),
		('coy', ctypes.POINTER(kinsol.realtype)),
		('ax', kinsol.realtype),
		('ay', kinsol.realtype),
		('dx', kinsol.realtype),
		('dy', kinsol.realtype),
		('uround', kinsol.realtype),
		('sqruround', kinsol.realtype),
		('mx', ctypes.c_long),
		('my', ctypes.c_long),
		('ns', ctypes.c_long),
		('np', ctypes.c_long)
	]
PUserData = ctypes.POINTER(UserData)

def WebRate(xx, yy, cxy, cxy_offset, ratesxy, ratesxy_offset, f_data):
	data = ctypes.cast(f_data, PUserData).contents
	
	for i in range(NUM_SPECIES):
		ratesxy[i+ratesxy_offset] = 0
		for j in range(NUM_SPECIES):
			ratesxy[i+ratesxy_offset] += cxy[cxy_offset+j] * data.acoef[i][j]
	
	fac = 1.0 + ALPHA * xx * yy
	
	for i in range(NUM_SPECIES):
		ratesxy[i+ratesxy_offset] = cxy[i+cxy_offset] * (data.bcoef[i] * fac + ratesxy[i+ratesxy_offset])

def func(cc, fval, f_data):
	data = ctypes.cast(f_data, PUserData).contents
	rates = kinsol.NVector(data.rates)
	
	for jy in range(MY):
		yy = data.dy*jy
		if jy != 0:
			idyl = NSMX
		else:
			idyl = -NSMX
		if jy != MY-1:
			idyu = NSMX
		else:
			idyu = -NSMX
		for jx in range(MX):
			xx = data.dx*jx
			if jx !=  0:
				idxl  = NUM_SPECIES 
			else:
				idxl  = -NUM_SPECIES
			if jx != MX-1:
				idxr  = NUM_SPECIES 
			else:
				idxr  = -NUM_SPECIES
			
			offset = jx*NUM_SPECIES+jy*NSMX
			
			WebRate(xx, yy, cc, offset, rates, offset, f_data)
			
			for i_s in range(NUM_SPECIES):
				dcyli = cc[offset+i_s] - cc[offset - idyl + i_s]
				dcyui = cc[offset + idyu + i_s] - cc[offset+i_s]
				dcxli = cc[offset+i_s] - cc[offset - idxl + i_s]
				dcxri = cc[offset + idxr + i_s] - cc[offset+i_s]
				fval[offset+i_s] = data.coy[i_s] * (dcyui - dcyli) + data.cox[i_s] * (dcxri - dcxli) + rates[offset+i_s]
	return 0

def PrecSetupBD(cc, cscale, fval, fscale, P_data, vtemp1, vtemp2):
	data = ctypes.cast(P_data, PUserData).contents
	rates = kinsol.NVector(data.rates)
	perturb_rates = kinsol.NVector([0]*NUM_SPECIES)
	
	fac = fval.wl2norm(fscale)
	r0 = THOUSAND * data.uround * fac * NEQ
	if r0 == 0.0:
		r0 = 1.0
	
	for jy in range(MY):
		yy = jy*data.dy
		for jx in range(MX):
			xx = jx*data.dx
			#Pxy = (data.P)[jx][jy]
			
			offset = jx*NUM_SPECIES+jy*NSMX

			for j in range(NUM_SPECIES): 
				csave = cc[offset+j]
				r = max([data.sqruround*abs(csave), r0/cscale[offset+j]])
				cc[offset+j] += r
				fac = 1.0/r
				
				WebRate(xx, yy, cc, offset, perturb_rates, 0, P_data)
				
				cc[offset+j] = csave
				
				for i in range(NUM_SPECIES):
					data.P[jx][jy][j][i] = (perturb_rates[i] - rates[offset+i]) * fac
			ret = kinsol.denGETRF(data.P[jx][jy], NUM_SPECIES, NUM_SPECIES, data.pivot[jx][jy])
			if ret != 0:
				return 1
	
	return 0

def PrecSolveBD(cc, cscale, fval, fscale, vv, P_data, ftem):
	data = ctypes.cast(P_data, PUserData).contents
	
	for jx in range(MX):
		for jy in range(MY):
			kinsol.denGETRS(data.P[jx][jy], NUM_SPECIES, data.pivot[jx][jy], vv.ptrto(jx*NUM_SPECIES + jy*NSMX))
	
	return 0

def PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol):
	print "\nPredator-prey test problem --  KINSol (serial version)\n"
	print "Mesh dimensions = %d X %d"%(MX, MY)
	print "Number of species = %d"%(NUM_SPECIES)
	print "Total system size = %d\n"%(NEQ)
	print "Flag globalstrategy = %d (0 = None, 1 = Linesearch)"%(globalstrategy)
	print "Linear solver is SPGMR with maxl = %d, maxlrst = %d"%(maxl, maxlrst)
	print "Preconditioning uses interaction-only block-diagonal matrix"
	print "Positivity constraints imposed on all components "
	print "Tolerance parameters:  fnormtol = %lg   scsteptol = %lg"%(fnormtol, scsteptol)
	
	print "\nInitial profile of concentration"
	print "At all mesh points:  %lg %lg %lg   %lg %lg %lg"%(PREYIN, PREYIN, PREYIN, PREDIN, PREDIN, PREDIN)
	

def PrintOutput(cc):
	jy = 0
	jx = 0
	offset = jx*NUM_SPECIES+jy*NSMX

	print "\nAt bottom left:",
	
	for i_s in range(NUM_SPECIES):
		if (i_s%6)*6 == i_s:
			print
		print"%lg "%(cc[offset+i_s])
	
	jy = MY-1
	jx = MX-1
	offset = jx*NUM_SPECIES+jy*NSMX

	print "\n\nAt top right:",
	
	for i_s in range(NUM_SPECIES):
		if (i_s%6)*6 == i_s:
			print
		print"%lg "%(cc[offset+i_s])
	
	print "\n"

def PrintFinalStats(kmem):
	nni = kinsol.KINGetNumNonlinSolvIters(kmem)
	nfe = kinsol.KINGetNumFuncEvals(kmem)
	nli = kinsol.KINSpilsGetNumLinIters(kmem)
	npe = kinsol.KINSpilsGetNumPrecEvals(kmem)
	nps = kinsol.KINSpilsGetNumPrecSolves(kmem)
	ncfl = kinsol.KINSpilsGetNumConvFails(kmem)
	nfeSG = kinsol.KINSpilsGetNumFuncEvals(kmem)
	
	print "Final Statistics.. "
	print "nni    = %5ld    nli   = %5ld"%(nni, nli)
	print "nfe    = %5ld    nfeSG = %5ld"%(nfe, nfeSG)
	print "nps    = %5ld    npe   = %5ld     ncfl  = %5ld"%(nps, npe, ncfl)

def SetInitialProfiles(cc, sc):
	ctemp = [PREYIN]*(NUM_SPECIES/2) + [PREDIN]*(NUM_SPECIES/2)
	stemp = [1]*(NUM_SPECIES/2) + [0.00001]*(NUM_SPECIES/2)
	
	for jy in range(MY):
		for jx in range(MX):
			offset = jx*NUM_SPECIES+jy*NSMX
			for i in range(NUM_SPECIES):
				cc[offset+i] = ctemp[i]
				sc[offset+i] = stemp[i]

globalstrategy = kinsol.KIN_NONE

data = UserData()
for jx in range(MX):
	for jy in range(MY):
		data.P[jx][jy] = kinsol.denalloc(NUM_SPECIES, NUM_SPECIES)
		data.pivot[jx][jy] = kinsol.denallocpiv(NUM_SPECIES)
data.acoef = kinsol.denalloc(NUM_SPECIES, NUM_SPECIES)
data.bcoef = (kinsol.realtype*NUM_SPECIES)()
data.cox = (kinsol.realtype*NUM_SPECIES)()
data.coy = (kinsol.realtype*NUM_SPECIES)()
rates = kinsol.NVector([0]*NEQ)
data.rates = rates.data

data.mx = MX
data.my = MY
data.ns = NUM_SPECIES
data.np = NUM_SPECIES/2
data.ax = AX
data.ay = AY
data.dx = (data.ax)/(MX-1)
data.dy = (data.ay)/(MY-1)
data.uround = nvecserial.UNIT_ROUNDOFF
data.sqruround = math.sqrt(data.uround)

dx2 = data.dx**2
dy2 = data.dy**2

for i in range(data.np):
    for j in range(data.np):
    	data.acoef[i][data.np+j] = -GG
    	data.acoef[i+data.np][0+j] = EE
    	data.acoef[i][0+j] = 0.0
    	data.acoef[i+data.np][data.np+j] = 0.0

    data.acoef[i][i] = -AA
    data.acoef[i+data.np][i+data.np] = -AA

    data.bcoef[i] = BB
    data.bcoef[i+data.np] = -BB

    data.cox[i] = DPREY/dx2
    data.cox[i+data.np] = DPRED/dx2

    data.coy[i] = DPREY/dy2
    data.coy[i+data.np] = DPRED/dy2

cc = kinsol.NVector([0]*NEQ)
sc = kinsol.NVector([0]*NEQ)
constraints = kinsol.NVector([2]*NEQ)

SetInitialProfiles(cc, sc)

kmem = kinsol.KINCreate()
kinsol.KINMalloc(kmem, func, cc)
kinsol.KINSetFdata(kmem, ctypes.pointer(data))
kinsol.KINSetConstraints(kmem, constraints)
kinsol.KINSetFuncNormTol(kmem, FTOL)
kinsol.KINSetScaledStepTol(kmem, STOL)

maxl = 15
maxlrst = 2
kinsol.KINSpgmr(kmem, maxl)
kinsol.KINSpilsSetMaxRestarts(kmem, maxlrst)
kinsol.KINSpilsSetPreconditioner(kmem, PrecSetupBD, PrecSolveBD, ctypes.pointer(data))

PrintHeader(globalstrategy, maxl, maxlrst, FTOL, STOL)

kinsol.KINSol(kmem, cc, globalstrategy, sc, sc)

print "\n\nComputed equilibrium species concentrations:"
PrintOutput(cc)

PrintFinalStats(kmem)

for jx in range(MX):
	for jy in range(MY):
		kinsol.denfree(data.P[jx][jy])
		kinsol.denfreepiv(data.pivot[jx][jy])

kinsol.denfree(data.acoef)
