#!/usr/bin/python

try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math
import sys

AA = 1.0
EE = 1.0e4
GG = 0.5e-6
BB = 1.0
DPREY = 1.0
DPRED = 0.5
ALPH = 1.0
NP = 3
NS = (2*NP)

MX = 6
MY = 6
MXNS = (MX*NS)
AX = 1.0
AY = 1.0
DX = (AX/float(MX-1))
DY = (AY/float(MY-1))
MP = NS
MQ = (MX*MY)
MXMP = (MX*MP)
NGX = 2
NGY = 2
NGRP = (NGX*NGY)
ITMAX = 5

NEQ = (NS*MX*MY)
T0 = 0.0
RTOL = 1.0e-5
ATOL = 1.0e-5

MAXL = 0
DELT = 0.0

T1 = 1.0e-8
TOUT_MULT = 10.0
DTOUT = 1.0
NOUT = 18

class WebData(ctypes.Structure):
	_fields_ = [
		("P", ctypes.POINTER(ctypes.POINTER(cvodes.realtype))*NGRP),
		("pivot", ctypes.POINTER(ctypes.c_long)*NGRP),
		("ns", ctypes.c_int),
		("mxns", ctypes.c_int),
		("mp", ctypes.c_int),
		("mq", ctypes.c_int),
		("mx", ctypes.c_int),
		("my", ctypes.c_int),
		("ngrp", ctypes.c_int),
		("ngx", ctypes.c_int),
		("ngy", ctypes.c_int),
		("mxmp", ctypes.c_int),
		("jgx", ctypes.c_int*(NGX+1)),
		("jgy", ctypes.c_int*(NGY+1)),
		("jigx", ctypes.c_int*(MX)),
		("jigy", ctypes.c_int*(MY)),
		("jxr", ctypes.c_int*(NGX)),
		("jyr", ctypes.c_int*(NGY)),
		("acoef", (cvodes.realtype*NS)*NS),
		("bcoef", cvodes.realtype*(NS)),
		("diff", cvodes.realtype*(NS)),
		("cox", cvodes.realtype*(NS)),
		("coy", cvodes.realtype*(NS)),
		("dx", cvodes.realtype),
		("dy", cvodes.realtype),
		("srur", cvodes.realtype),
		("fsave", cvodes.realtype*(NEQ)),
		("rewt", nvecserial.PVector),
		("cvodes_mem", ctypes.c_void_p)
	]
PWebData = ctypes.POINTER(WebData)

def InitUserData(wdata):
	for j in range(NS):
		for i in range(NS):
			wdata.acoef[i][j] = 0.0
	
	for j in range(NP):
		for i in range(NP):
			wdata.acoef[NP+i][j] = EE
			wdata.acoef[i][NP+j] = -GG
		wdata.acoef[j][j] = -AA
		wdata.acoef[NP+j][NP+j] = -AA
		wdata.bcoef[j] = BB
		wdata.bcoef[NP+j] = -BB
		wdata.diff[j] = DPREY
		wdata.diff[NP+j] = DPRED
	
	wdata.ns = NS
	wdata.mxns = MXNS
	dx = wdata.dx = DX
	dy = wdata.dy = DY
	for i in range(NS):
		wdata.cox[i] = wdata.diff[i]/(dx**2)
		wdata.coy[i] = wdata.diff[i]/(dy**2)
	
	wdata.mp = MP
	wdata.mq = MQ
	wdata.mx = MX
	wdata.my = MY
	wdata.srur = math.sqrt(nvecserial.UNIT_ROUNDOFF)
	wdata.mxmp = MXMP
	wdata.ngrp = NGRP
	wdata.ngx = NGX
	wdata.ngy = NGY
	SetGroups(MX, NGX, wdata.jgx, wdata.jigx, wdata.jxr)
	SetGroups(MY, NGY, wdata.jgy, wdata.jigy, wdata.jyr)

def SetGroups(m, ng, jg, jig, jr):
	mper = m/ng;
	for ig in range(ng):
		jg[ig] = ig*mper
	jg[ng] = m
	
	ngm1 = ng - 1
	len1 = ngm1*mper
	for j in range(len1):
		jig[j] = j/mper
	for j in range(len1, m):
		jig[j] = ngm1
	
	for ig in range(ngm1):
		jr[ig] = ((2*ig+1)*mper-1)/2
	jr[ngm1] = (ngm1*mper+m-1)/2

def CInit(c, wdata):
	ns = wdata.ns
	mxns = wdata.mxns
	dx = wdata.dx
	dy = wdata.dy
	
	x_factor = 4.0/(AX**2)
	y_factor = 4.0/(AY**2)
	for jy in range(MY):
		y = jy*dy
		argy = (y_factor*y*(AY-y))**2
		iyoff = mxns*jy
		for jx in range(MX): 
			x = jx*dx
			argx = (x_factor*x*(AX-x))**2
			ioff = iyoff + ns*jx
			for i in range(1,ns+1):
				ici = ioff + i-1
				c[ici] = 10.0 + i*argx*argy

def PrintIntro():
	print "\n\nDemonstration program for CVODE - CVSPGMR linear solver\n"
	print "Food web problem with ns species, ns = %d"%(NS)
	print "Predator-prey interaction and diffusion on a 2-D square\n"
	print "Matrix parameters: a = %.2lg   e = %.2lg   g = %.2lg"%(AA, EE, GG)
	print "b parameter = %.2lg"%(BB)
	print "Diffusion coefficients: Dprey = %.2lg   Dpred = %.2lg"%(DPREY, DPRED)
	print "Rate parameter alpha = %.2lg\n"%(ALPH)
	print "Mesh dimensions (mx,my) are %d, %d. "%(MX, MY),
	print "Total system size is neq = %d \n"%(NEQ)
	print "Tolerances: itol = %s,  reltol = %.2lg, abstol = %.2lg \n"%("CV_SS", RTOL, ATOL)
	print "Preconditioning uses a product of:"
	print "  (1) Gauss-Seidel iterations with",
	print "itmax = %d iterations, and"%(ITMAX)
	print "  (2) interaction-only block-diagonal matrix",
	print "with block-grouping"
	print "  Number of diagonal block groups = ngrp = %d"%(NGRP),
	print " (ngx by ngy, ngx = %d, ngy = %d)"%(NGX, NGY)
	print "\n\n----------------------------------------------------------------------------"

def PrintHeader(jpre, gstype):
	if jpre == cvodes.PREC_LEFT:
		print "\n\nPreconditioner type is           jpre = %s"%("PREC_LEFT")
	else:
		print "\n\nPreconditioner type is           jpre = %s"%("PREC_RIGHT")
	
	if	gstype == cvodes.MODIFIED_GS:
		print "\nGram-Schmidt method type is    gstype = %s\n\n"%("MODIFIED_GS")
	else:
		print "\nGram-Schmidt method type is    gstype = %s\n\n"%("CLASSICAL_GS")

def PrintAllSpecies(c, ns, mxns, t):
	print "c values at t = %lg:\n"%(t.value)
	for i in range(1, ns+1):
		print "Species %d"%(i)
		for jy in range(MY-1,-1,-1):
			for jx in range(MX):
				sys.stdout.write("%-10.6lg"%(c[(i-1) + jx*ns + jy*mxns]))
			print
		print

def PrintOutput(cvodes_mem, t):
	nst = cvodes.CVodeGetNumSteps(cvodes_mem)
	nfe = cvodes.CVodeGetNumRhsEvals(cvodes_mem)
	nni = cvodes.CVodeGetNumNonlinSolvIters(cvodes_mem)
	qu = cvodes.CVodeGetLastOrder(cvodes_mem)
	hu = cvodes.CVodeGetLastStep(cvodes_mem)
	
	print "t = %10.2le  nst = %ld  nfe = %ld  nni = %ld"%(t.value, nst, nfe, nni),
	print " qu = %d  hu = %11.2le\n"%(qu, hu)

def PrintFinalStats(cvodes_mem):
	(lenrw, leniw) = cvodes.CVodeGetWorkSpace(cvodes_mem)
	nst = cvodes.CVodeGetNumSteps(cvodes_mem)
	nfe = cvodes.CVodeGetNumRhsEvals(cvodes_mem)
	nsetups = cvodes.CVodeGetNumLinSolvSetups(cvodes_mem)
	netf = cvodes.CVodeGetNumErrTestFails(cvodes_mem)
	nni = cvodes.CVodeGetNumNonlinSolvIters(cvodes_mem)
	ncfn = cvodes.CVodeGetNumNonlinSolvConvFails(cvodes_mem)
	
	(lenrwLS, leniwLS) = cvodes.CVSpilsGetWorkSpace(cvodes_mem)
	nli = cvodes.CVSpilsGetNumLinIters(cvodes_mem)
	npe = cvodes.CVSpilsGetNumPrecEvals(cvodes_mem)
	nps = cvodes.CVSpilsGetNumPrecSolves(cvodes_mem)
	ncfl = cvodes.CVSpilsGetNumConvFails(cvodes_mem)
	nfeLS = cvodes.CVSpilsGetNumRhsEvals(cvodes_mem)
	
	print "\n\n Final statistics for this run:\n"
	print " CVode real workspace length           = %4ld "%(lenrw.value)
	print " CVode integer workspace length        = %4ld "%(leniw.value)
	print " CVSPGMR real workspace length         = %4ld "%(lenrwLS)
	print " CVSPGMR integer workspace length      = %4ld "%(leniwLS)
	print " Number of steps                       = %4ld "%(nst)
	print " Number of f-s                         = %4ld "%(nfe)
	print " Number of f-s (SPGMR)                 = %4ld "%(nfeLS)
	print " Number of f-s (TOTAL)                 = %4ld "%(nfe + nfeLS)
	print " Number of setups                      = %4ld "%(nsetups)
	print " Number of nonlinear iterations        = %4ld "%(nni)
	print " Number of linear iterations           = %4ld "%(nli)
	print " Number of preconditioner evaluations  = %4ld "%(npe)
	print " Number of preconditioner solves       = %4ld "%(nps)
	print " Number of error test failures         = %4ld "%(netf)
	print " Number of nonlinear conv. failures    = %4ld "%(ncfn)
	print " Number of linear convergence failures = %4ld "%(ncfl)
	if nni > 0:
		avdim = float(nli)/nni
	else:
		avdim = 0.0
	print " Average Krylov subspace dimension     = %.3f "%(avdim)
	print "\n\n----------------------------------------------------------------------------"
	print     "----------------------------------------------------------------------------"

def WebRates(x, y, t, c, c_off, rate, rate_off, wdata):
	ns = wdata.ns
	acoef = wdata.acoef
	bcoef = wdata.bcoef
	
	for i in range(ns):
		rate[i+rate_off] = 0.0
	
	for j in range(ns):
		for i in range(ns):
			rate[i+rate_off] += c[j+c_off] * acoef[i][j]
	
	fac = 1.0 + ALPH*x*y
	for i in range(ns):
		rate[i+rate_off] = c[i+c_off]*(bcoef[i]*fac + rate[i+rate_off])

def f(t, c, cdot, f_data):
	wdata = ctypes.cast(f_data, PWebData).contents
	
	mxns = wdata.mxns
	ns = wdata.ns
	fsave = wdata.fsave
	cox = wdata.cox
	coy = wdata.coy
	mxns = wdata.mxns
	dx = wdata.dx
	dy = wdata.dy
	
	for jy in range(MY):
		y = jy*dy
		iyoff = mxns*jy
		if jy == MY-1:
			idyu = -mxns
		else:
			idyu = mxns
		if jy == 0:
			idyl = -mxns
		else:
			idyl = mxns
		for jx in range(MX):
			x = jx*dx
			ic = iyoff + ns*jx
			WebRates(x, y, t, c, ic, fsave, ic, wdata)
			if jx == MX-1:
				idxu = -ns
			else:
				idxu = ns
			if jx == 0:
				idxl = -ns
			else:
				idxl = ns
			for i in range(1, ns+1):
				ici = ic + i-1
				dcyli = c[ici] - c[ici-idyl]
				dcyui = c[ici+idyu] - c[ici]
				dcxli = c[ici] - c[ici-idxl]
				dcxui = c[ici+idxu] - c[ici]
				cdot[ici] = coy[i-1]*(dcyui - dcyli) + cox[i-1]*(dcxui - dcxli) + fsave[ici]
	
	return 0

def fblock(t, c, jx, jy, cdot, wdata):
	iblok = jx + jy*(wdata.mx)
	y = jy*(wdata.dy)
	x = jx*(wdata.dx)
	ic = (wdata.ns)*(iblok)
	WebRates(x, y, t, c, ic, cdot, 0, wdata)

def Precond(t, c, fc, jok, jcurPtr, gamma, P_data, vtemp1, vtemp2, vtemp3):
	wdata = ctypes.cast(P_data, PWebData).contents
	cvodes_mem = wdata.cvodes_mem
	rewt = nvecserial.NVector(wdata.rewt)
	cvodes.CVodeGetErrWeights(cvodes_mem, rewt)
	
	uround = nvecserial.UNIT_ROUNDOFF
	
	P = wdata.P
	pivot = wdata.pivot
	jxr = wdata.jxr
	jyr = wdata.jyr
	mp = wdata.mp
	srur = wdata.srur
	ngrp = wdata.ngrp
	ngx = wdata.ngx
	ngy = wdata.ngy
	mxmp = wdata.mxmp
	fsave = wdata.fsave
	
	fac = fc.wrmsnorm(rewt)
	r0 = 1000.0*abs(gamma)*uround*NEQ*fac
	if r0 == 0.0:
		r0 = 1.0
	
	for igy in range(ngy):
		jy = jyr[igy]
		if00 = jy*mxmp
		for igx in range(ngx):
			jx = jxr[igx]
			if0 = if00 + jx*mp
			ig = igx + igy*ngx
			for j in range(mp):
				jj = if0 + j
				save = c[jj]
				r = max([srur*abs(save),r0/rewt[jj]])
				c[jj] += r
				fac = -gamma/r
				fblock (t, c, jx, jy, vtemp1, wdata)
				for i in range(mp):
					P[ig][j][i] = (vtemp1[i] - fsave[if0+i])*fac
				c[jj] = save
	
	for ig in range(ngrp):
		cvodes.denaddI(P[ig], mp)
		ier = cvodes.denGETRF(P[ig], mp, mp, pivot[ig])
		if ier != 0:
			return 1
	
	jcurPtr.contents.value = 1
	return 0

def v_inc_by_prod(u, u_off, v, v_off, w, w_off, n):
	for i in range(n):
		u[i+u_off] += v[i+v_off]*w[i+w_off]

def v_sum_prods(u, u_off, p, p_off, q, q_off, v, v_off, w, w_off, n):
	for i in range(n):
		u[i+u_off] = p[i+p_off]*q[i+q_off] + v[i+v_off]*w[i+w_off]

def v_prod(u, u_off, v, v_off, w, w_off, n):
	for i in range(n):
		u[i+u_off] = v[i+v_off]*w[i+w_off]

def v_zero(u, u_off, n):
	for i in range(n):
		u[i+u_off] = 0.0

def GSIter(gamma, z, x, wdata):
	beta = [0]*NS
	beta2 = [0]*NS
	cof1 = [0]*NS
	gam = [0]*NS
	gam2 = [0]*NS
	
	ns = wdata.ns
	mx = wdata.mx
	my = wdata.my
	mxns = wdata.mxns
	cox = wdata.cox
	coy = wdata.coy
	
	for i in range(ns):
		temp = 1.0/(1.0 + 2.0*gamma*(cox[i] + coy[i]))
		beta[i] = gamma*cox[i]*temp
		beta2[i] = 2.0*beta[i]
		gam[i] = gamma*coy[i]*temp
		gam2[i] = 2.0*gam[i]
		cof1[i] = temp
	
	for jy in range(my):
		iyoff = mxns*jy
		for jx in range(mx):
			ic = iyoff + ns*jx
			v_prod(x, ic, cof1, 0, z, ic, ns)
	z[:] = [0]*len(z)
	
	for iter in range(1, ITMAX+1):
		if iter > 1:
			for jy in range(my):
				iyoff = mxns*jy
				for jx in range(mx):
					ic = iyoff + ns*jx
					if jx == 0:
						x_loc = 0
					else:
						if jx == mx-1:
							x_loc = 2
						else:
							xloc = 1
					if jy == 0:
						y_loc = 0
					else:
						if jy == my-1:
							y_loc = 2
						else:
							yloc = 1
					if (3*y_loc+x_loc) == 0:
						v_sum_prods(x, ic, beta2, 0, x, ic+ns, gam2, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 1: 
						v_sum_prods(x, ic, beta, 0, x, ic+ns, gam2, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 2: 
						v_prod(x, ic, gam2, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 3: 
						v_sum_prods(x, ic, beta2, 0, x, ic+ns, gam, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 4: 
						v_sum_prods(x, ic, beta, 0, x, ic+ns, gam, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 5: 
						v_prod(x, ic, gam, 0, x, ic+mxns, ns)
					elif (3*y_loc+x_loc) == 6: 
						v_prod(x, ic, beta2, 0, x, ic+ns, ns)
					elif (3*y_loc+x_loc) == 7: 
						v_prod(x, ic, beta, 0, x, ic+ns, ns)
					elif (3*y_loc+x_loc) == 8: 
						v_zero(x, ic, ns)
		
		for jy in range(my):
			iyoff = mxns*jy
			for jx in range(mx):
				ic = iyoff + ns*jx
				if jx == 0:
					x_loc = 0
				else:
					if jx == mx-1:
						x_loc = 2
					else:
						xloc = 1
				if jy == 0:
					y_loc = 0
				else:
					if jy == my-1:
						y_loc = 2
					else:
						yloc = 1
				if (3*y_loc+x_loc) == 0: 
				  pass
				elif (3*y_loc+x_loc) == 1: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-ns, ns)
				elif (3*y_loc+x_loc) == 2: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-ns, ns)
				elif (3*y_loc+x_loc) == 3: 
				  v_inc_by_prod(x, ic, gam, 0, x, ic-mxns, ns)
				elif (3*y_loc+x_loc) == 4: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-ns, ns)
				  v_inc_by_prod(x, ic, gam, 0, x, ic-mxns, ns)
				elif (3*y_loc+x_loc) == 5: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-ns, ns)
				  v_inc_by_prod(x, ic, gam, 0, x, ic-mxns, ns)
				elif (3*y_loc+x_loc) == 6: 
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-mxns, ns)
				elif (3*y_loc+x_loc) == 7: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-ns, ns)
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-mxns, ns)
				elif (3*y_loc+x_loc) == 8: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-ns, ns)
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-mxns, ns)
		
		z[:] = z.linearsum(1.0, 1.0, x)

def PSolve(tn, c, fc, r, z, gamma, delta, lr, P_data, vtemp):
	wdata = ctypes.cast(P_data, PWebData).contents
	
	z[:] = r
	
	GSIter(gamma, z, vtemp, wdata)
	
	P = wdata.P
	pivot = wdata.pivot
	mx = wdata.mx
	my = wdata.my
	ngx = wdata.ngx
	mp = wdata.mp
	jigx = wdata.jigx
	jigy = wdata.jigy
	
	iv = 0
	for jy in range(my):
		igy = jigy[jy]
		for jx in range(mx):
			igx = jigx[jx]
			ig = igx + igy*ngx
			cvodes.denGETRS(P[ig], mp, pivot[ig], z.ptrto(iv))
			iv += mp
	
	return 0

c = cvodes.NVector([0]*NEQ)
rewt = cvodes.NVector([0]*NEQ)
wdata = WebData()
for i in range(NGRP):
	wdata.P[i] = cvodes.denalloc(NS,NS)
	wdata.pivot[i] = cvodes.denallocpiv(NS)
wdata.rewt = rewt.data
InitUserData(wdata)
ns = wdata.ns
mxns = wdata.mxns

PrintIntro()

for jpre in range(cvodes.PREC_LEFT, cvodes.PREC_RIGHT+1):
	for gstype in range(cvodes.MODIFIED_GS,cvodes.CLASSICAL_GS+1):
		t = cvodes.realtype(T0)
		CInit(c, wdata)
		PrintHeader(jpre, gstype)
		
		firstrun = (jpre == cvodes.PREC_LEFT and gstype == cvodes.MODIFIED_GS)
		if firstrun:
			cvodes_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
			wdata.cvodes_mem = cvodes_mem.obj
			cvodes.CVodeSetFdata(cvodes_mem, ctypes.pointer(wdata))
			cvodes.CVodeMalloc(cvodes_mem, f, t.value, c, cvodes.CV_SS, RTOL, ATOL)
			cvodes.CVSpgmr(cvodes_mem, jpre, MAXL)
			cvodes.CVSpilsSetGSType(cvodes_mem, gstype)
			cvodes.CVSpilsSetDelt(cvodes_mem, DELT)
			cvodes.CVSpilsSetPreconditioner(cvodes_mem, Precond, PSolve, ctypes.pointer(wdata))
		
			PrintAllSpecies(c, ns, mxns, t)
		else:
			cvodes.CVodeReInit(cvodes_mem, f, t.value, c, cvodes.CV_SS, RTOL, ATOL)
			cvodes.CVSpilsSetPrecType(cvodes_mem, jpre)
			cvodes.CVSpilsSetGSType(cvodes_mem, gstype)
		
		tout = T1
		iout = 1
		while iout <= NOUT:
			cvodes.CVode(cvodes_mem, tout, c, ctypes.byref(t), cvodes.CV_NORMAL)
			PrintOutput(cvodes_mem, t)
			sys.stderr.write("t = %g\n"%(t.value))
			if firstrun and iout % 3 == 0:
				PrintAllSpecies(c, ns, mxns, t)
			if tout > 0.9:
				tout += DTOUT
			else:
				tout *= TOUT_MULT
			iout += 1
		
		PrintFinalStats(cvodes_mem)

for i in range(wdata.ngrp):
	cvodes.denfree(wdata.P[i])
	cvodes.denfreepiv(wdata.pivot[i])
