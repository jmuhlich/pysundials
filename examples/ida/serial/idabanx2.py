try:
	from pysundials import ida
	from pysundials import nvecserial
except ImportError:
	import ida
	import nvecserial
import math
import ctypes

NPREY = 1
NUM_SPECIES = 2*NPREY

FOURPI = 4.0*math.pi

MX = 20
MY = 20
NSMX = (NUM_SPECIES * MX)
NEQ = (NUM_SPECIES*MX*MY)
AA = 1.0
EE = 10000.
GG = 0.5e-6
BB = 1.0
DPREY = 1.0
DPRED = 0.05
ALPHA = 50.
BETA = 1000.
AX = 1.0
AY = 1.0
RTOL = 1.e-5
ATOL = 1.e-5
NOUT = 6
TMULT = 10.0
TADD = 0.3

#define IJ_Vptr(vv,i,j) (&NV_Ith_S(vv, (i)*NUM_SPECIES + (j)*NSMX))

class UserData(ctypes.Structure):
	_fields_ = [
		('Neq', ctypes.c_long),
		('ns', ctypes.c_long),
		('np', ctypes.c_long),
		('mx', ctypes.c_long),
		('my', ctypes.c_long),
		('dx', ida.realtype), 
		('dy', ida.realtype), 
		('acoef', ctypes.POINTER(ctypes.POINTER(ida.realtype))),
		('cox', ida.realtype*NUM_SPECIES),
		('coy', ida.realtype*NUM_SPECIES),
		('bcoef', ida.realtype*NUM_SPECIES),
		('rates', nvecserial.PVector)
	]
PUserData = ctypes.POINTER(UserData)

def WebRates(xx, yy, cc, webdata, offset):
	r = ida.NVector(webdata.rates)
	for i_s in range(NUM_SPECIES):
		r[offset+i_s] = 0
		for i in range(NUM_SPECIES):
			r[offset+i_s] += cc[offset+i] * webdata.acoef[i_s][i]
	
	fac = 1.0 + ALPHA*xx*yy + BETA*math.sin(FOURPI*xx)*math.sin(FOURPI*yy)
	
	for i_s in range(NUM_SPECIES):
		r[offset+i_s] = cc[offset+i_s]*(webdata.bcoef[i_s]*fac + r[offset+i_s])

def Fweb(tcalc, cc, crate, webdata):
	r = ida.NVector(webdata.rates)
	for jy in range(MY):
		yy = (webdata.dy) * jy 
		if jy != MY-1:
			idyu = NSMX
		else:
			idyu = -NSMX
		if jy != 0:
			idyl = NSMX
		else:
			idyl = -NSMX
		
		for jx in range(MX):
			xx = (webdata.dx) * jx
			if jx != MX-1:
				idxu = NUM_SPECIES
			else:
				idxu = -NUM_SPECIES
			if jx != 0: 
				idxl = NUM_SPECIES
			else:
				idxl = -NUM_SPECIES
			offset = jx*NUM_SPECIES + jy*NSMX
			
			WebRates(xx, yy, cc, webdata, offset)
			
			for i_s in range(NUM_SPECIES):
				dcyli = cc[offset+i_s] - cc[offset - idyl + i_s]
				dcyui = cc[offset + idyu + i_s] - cc[offset+i_s]
				dcxli = cc[offset+i_s] - cc[offset - idxl + i_s]
				dcxui = cc[offset + idxu + i_s] - cc[offset+i_s]
				crate[offset+i_s] = webdata.coy[i_s] * (dcyui - dcyli) + webdata.cox[i_s] * (dcxui - dcxli) + r[offset+i_s]
				
def resweb(tt, cc, cp, res,	rdata):
	webdata = ctypes.cast(rdata, PUserData).contents
	
	Fweb(tt, cc, res, webdata)
	
	for jy in range(MY):
		yloc = NSMX * jy
		for jx in range(MX):
			loc = yloc + NUM_SPECIES * jx
			for i_s in range(NUM_SPECIES):
				if (i_s < webdata.np):
					res[loc+i_s] = cp[loc+i_s] - res[loc+i_s]
				else:
					res[loc+i_s] = -res[loc+i_s]
	return 0

def InitUserData(webdata):
	webdata.mx = MX
	webdata.my = MY
	webdata.ns = NUM_SPECIES
	webdata.np = NPREY
	webdata.dx = AX/(MX-1)
	webdata.dy = AY/(MY-1)
	webdata.Neq= NEQ
	
	for i in range(webdata.np):
		for j in range(webdata.np):
			webdata.acoef[i][webdata.np+j] = -GG
			webdata.acoef[i+webdata.np][0+j] = EE
			webdata.acoef[i][0+j]= 0.0
			webdata.acoef[i+webdata.np][webdata.np+j]= 0.0
		
		webdata.acoef[i][i]= -AA
		webdata.acoef[i+webdata.np][i+webdata.np]= -AA
		
		webdata.bcoef[i] = BB
		webdata.bcoef[i+webdata.np] = -BB
		webdata.cox[i] = DPREY/(webdata.dx**2)
		webdata.cox[i+webdata.np] = DPRED/(webdata.dx**2)
		webdata.coy[i] = DPREY/(webdata.dy**2)
		webdata.coy[i+webdata.np] = DPRED/(webdata.dy**2)

def SetInitialProfiles(cc, cp, id, webdata):
	for jy in range(MY):
		yy = jy * webdata.dy
		yloc = NSMX * jy
		for jx in range(MX):
			xx = jx * webdata.dx
			xyfactor = 16.0*xx*(1.0-xx)*yy*(1.0-yy)
			xyfactor *= xyfactor
			loc = yloc + NUM_SPECIES*jx
			fac = 1.0 + ALPHA * xx * yy + BETA * math.sin(FOURPI*xx) * math.sin(FOURPI*yy)
			
			for i_s in range(NUM_SPECIES):
				if i_s < webdata.np:
					cc[loc+i_s] = 10.0 + (i_s+1) * xyfactor
					id[loc+i_s] = 1.0
				else:
					cc[loc+i_s] = 1.0e5
					id[loc+i_s] = 0.0
	
	Fweb(0.0, cc, cp, webdata)
	
	for jy in range(MY):
		yloc = NSMX * jy
		for jx in range(MX):
			loc = yloc + NUM_SPECIES * jx
			for i_s in range(webdata.np, NUM_SPECIES):
				cp[loc+i_s] = 0.0

def PrintHeader(mu, ml, rtol, atol):
	print "\nidabanx2: Predator-prey DAE serial example problem for IDA \n"
	print "Number of species ns: %d"%(NUM_SPECIES),
	print "		 Mesh dimensions: %d x %d"%(MX, MY),
	print "		 System size: %d"%(NEQ)
	print "Tolerance parameters:	rtol = %lg	 atol = %lg"%(rtol, atol)
	print "Linear solver: IDABAND,	Band parameters mu = %ld, ml = %ld"%(mu, ml)
	print "CalcIC called to correct initial predator concentrations.\n"
	print "-----------------------------------------------------------"
	print "	t				bottom-left	top-right 		| nst	k			h"
	print "-----------------------------------------------------------\n"

def PrintOutput(mem, c, t):
	kused = ida.IDAGetLastOrder(mem)
	nst = ida.IDAGetNumSteps(mem)
	hused = ida.IDAGetLastStep(mem)
	
	print "%8.2le %12.4le %12.4le	 | %3ld	%1d %12.4le"%(t, c[0], c[((MX-1)*NUM_SPECIES + (MY-1)*NSMX)+1], nst, kused, hused)
	for i in range(1, NUM_SPECIES):
		print "				 %12.4le %12.4le	 |"%(c[i], c[((MX-1)*NUM_SPECIES + (MY-1)*NSMX)+i])

def PrintFinalStats(mem):
	nst = ida.IDAGetNumSteps(mem)
	nni = ida.IDAGetNumNonlinSolvIters(mem)
	nre = ida.IDAGetNumResEvals(mem)
	netf = ida.IDAGetNumErrTestFails(mem)
	ncfn = ida.IDAGetNumNonlinSolvConvFails(mem)
	nje = ida.IDABandGetNumJacEvals(mem)
	nreLS = ida.IDABandGetNumResEvals(mem)

	print "-----------------------------------------------------------"
	print "Final run statistics: \n"
	print "Number of steps										= %ld"%(nst)
	print "Number of residual evaluations		 = %ld"%(nre+nreLS)
	print "Number of Jacobian evaluations		 = %ld"%(nje)
	print "Number of nonlinear iterations		 = %ld"%(nni)
	print "Number of error test failures			= %ld"%(netf)
	print "Number of nonlinear conv. failures = %ld"%(ncfn)

wrates = ida.NVector([0]*NEQ)
webdata = UserData()
webdata.rates = wrates.data
webdata.acoef = ida.denalloc(NUM_SPECIES, NUM_SPECIES)

InitUserData(webdata)

cc = ida.NVector([0]*NEQ)
cp = ida.NVector([0]*NEQ)
id = ida.NVector([0]*NEQ)

SetInitialProfiles(cc, cp, id, webdata)

t0 = 0.0

mem = ida.IDACreate()
ida.IDASetRdata(mem, ctypes.pointer(webdata))
ida.IDASetId(mem, id)
ida.IDAMalloc(mem, resweb, t0, cc, cp, ida.IDA_SS, RTOL, ATOL)
ida.IDABand(mem, NEQ, NSMX, NSMX)

tout = 0.001
ida.IDACalcIC(mem, ida.IDA_YA_YDP_INIT, tout)

PrintHeader(NSMX, NSMX, RTOL, ATOL)
PrintOutput(mem, cc, 0.0)

tret = ida.realtype(0.0)
for iout in range(1,NOUT+1):
	ida.IDASolve(mem, tout, ctypes.byref(tret), cc, cp, ida.IDA_NORMAL)
	PrintOutput(mem, cc, tret.value)
	if iout < 3:
		tout *= TMULT
	else:
		tout += TADD

PrintFinalStats(mem)
