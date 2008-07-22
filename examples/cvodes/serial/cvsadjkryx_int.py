#!/usr/bin/python

try:
	from pysundials import cvodes
	from pysundials import nvecserial
except ImportError:
	import cvodes
	import nvecserial
import ctypes
import math

AA = 1.0
EE = 1.0e4
GG = 0.5e-6
BB = 1.0
DPREY = 1.0
DPRED = 0.5
ALPH = 1.0
NP = 3
NS = (2*NP)

MX = 20
MY = 20
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

TOUT = 10.0
NSTEPS = 300
ISPEC = 6

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
		("fBsave", cvodes.realtype*(NEQ)),
		("rewt", nvecserial.PVector),
		("cvadj_mem", ctypes.c_void_p),
		("cvode_memF", ctypes.c_void_p)
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
	mper = m/ng
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
	c[NEQ] = 0.0

def PrintOutput(cB, ns, mxns, wdata):
	x = y = 0.0
	
	for i in range (ns+1):
		cmax = 0.0
		
		for jy in range(MY-1, -1, -1):
			for jx in range(MX):
				cij = c[(i-1) + jx*ns + jy*mxns]
				if abs(cij) > cmax:
					cmax = cij
					x = jx * wdata.dx
					y = jy * wdata.dy
		
		print "\nMaximum sensitivity with respect to I.C. of species %d"%(i)
		print "  lambda max = %le"%(cmax)
		print "at"
		print "  x = %le\n  y = %le"%(x, y)

def doubleIntgr(c, i, wdata):
	jy = 0
	intgr_x = c[(i-1)+jy*wdata.mxns]
	for jx in range(1,wdata.mx):
		intgr_x += 2.0*c[(i-1) + jx*wdata.ns + jy*wdata.mxns] 
	intgr_x += c[(i-1)+(wdata.mx-1)*wdata.ns+jy*wdata.mxns]
	intgr_x *= 0.5*wdata.dx
	
	intgr_xy = intgr_x
	
	for jy in range(1,wdata.my-1):
		intgr_x = c[(i-1)+jy*wdata.mxns]
		for jx in range(1,wdata.mx-1):
			intgr_x += 2.0*c[(i-1) + jx*wdata.ns + jy*wdata.mxns]
		intgr_x += c[(i-1)+(wdata.mx-1)*wdata.ns+jy*wdata.mxns]
		intgr_x *= 0.5*wdata.dx
		intgr_xy += 2.0*intgr_x
	
	jy = wdata.my-1
	intgr_x = c[(i-1)+jy*wdata.mxns]
	for jx in range(1,wdata.mx-1):
		intgr_x += 2.0*c[(i-1) + jx*wdata.ns + jy*wdata.mxns]
		
	intgr_x += c[(i-1)+(wdata.mx-1)*wdata.ns+jy*wdata.mxns]
	intgr_x *= 0.5*wdata.dx
	intgr_xy += intgr_x
	intgr_xy *= 0.5*wdata.dy

	return intgr_xy

def WebRates(x, y, t, c, c_off, rate, rate_off, wdata):
	for i in range(wdata.ns):
		rate[i+rate_off] = 0.0
	
	for j in range(wdata.ns):
		for i in range(wdata.ns):
			rate[i+rate_off] += c[j+c_off] * wdata.acoef[i][j]
	
	fac = 1.0 + ALPH*x*y
	for i in range(wdata.ns):
		rate[i+rate_off] = c[i+c_off]*(wdata.bcoef[i]*fac + rate[i+rate_off])

def WebRatesB(x, y, t, c, c_off, cB, cB_off, rate, rate_off, rateB, rateB_off, wdata):
	fac = 1.0 + ALPH*x*y
	
	for i in range(wdata.ns):
		rate[i+rate_off] = wdata.bcoef[i]*fac
	
	for j in range(wdata.ns):
		for i in range(wdata.ns):
			rate[i+rate_off] += wdata.acoef[i][j]*c[j+c_off]
	
	for i in range(wdata.ns):
		rateB[i+rateB_off] = cB[i+cB_off]*rate[i+rate_off]
		rate[i+rate_off] = c[i+c_off]*rate[i+rate_off]
	
	for j in range(wdata.ns):
		for i in range(wdata.ns):
			rateB[i+rateB_off] += wdata.acoef[j][i]*c[j+c_off]*cB[j+cB_off]

def f(t, c, cdot, f_data):
	wdata = ctypes.cast(f_data, PWebData).contents
	
	for jy in range(MY):
		y = jy*wdata.dy
		iyoff = wdata.mxns*jy
		if jy == MY-1:
			idyu = -wdata.mxns
		else:
			idyu = wdata.mxns
		if jy == 0:
			idyl = -wdata.mxns
		else:
			idyl = wdata.mxns
		for jx in range(MX):
			x = jx*wdata.dx
			ic = iyoff + wdata.ns*jx
			WebRates(x, y, t, c, ic, wdata.fsave, ic, wdata)
			if jx == MX-1:
				idxu = -wdata.ns
			else:
				idxu = wdata.ns
			if jx == 0:
				idxl = -wdata.ns
			else:
				idxl = wdata.ns
			for i in range(1, wdata.ns+1):
				ici = ic + i-1
				dcyli = c[ici] - c[ici-idyl]
				dcyui = c[ici+idyu] - c[ici]
				dcxli = c[ici] - c[ici-idxl]
				dcxui = c[ici+idxu] - c[ici]
				cdot[ici] = wdata.coy[i-1]*(dcyui - dcyli) + wdata.cox[i-1]*(dcxui - dcxli) + wdata.fsave[ici]
		cdot[NEQ] = doubleIntgr(c,ISPEC,wdata)
	
	return 0

def fB(t, c, cB, cBdot, f_data):
	wdata = ctypes.cast(f_data, PWebData).contents
	
	gu = [0]*NS
	gu[ISPEC-1] = 1.0
	
	for jy in range(MY):
		y = jy*wdata.dy
		iyoff = wdata.mxns*jy
		if jy == MY-1:
			idyu = -wdata.mxns
		else:
			idyu = wdata.mxns
		if jy == 0:
			idyl = -wdata.mxns
		else:
			idyl = wdata.mxns
		for jx in range(MX):
			x = jx*wdata.dx
			ic = iyoff + wdata.ns*jx
			WebRates(x, y, t, c, ic, cB, ic, wdata.fsave, ic, wdata.fBsave, ic, wdata)
			if jx == MX-1:
				idxu = -wdata.ns
			else:
				idxu = wdata.ns
			if jx == 0:
				idxl = -wdata.ns
			else:
				idxl = wdata.ns
			for i in range(1, wdata.ns+1):
				ici = ic + i-1
				dcyli = cB[ici] - cB[ici-idyl]
				dcyui = cB[ici+idyu] - cB[ici]
				dcxli = cB[ici] - cB[ici-idxl]
				dcxui = cB[ici+idxu] - cB[ici]
				cBdot[ici] = - wdata.coy[i-1]*(dcyui - dcyli) - wdata.cox[i-1]*(dcxui - dcxli) - wdata.fBsave[ici] - gu[i-1]

	return 0

def fblock(t, c, jx, jy, cdot, wdata):
	iblok = jx + jy*(wdata.mx)
	y = jy*(wdata.dy)
	x = jx*(wdata.dx)
	ic = (wdata.ns)*(iblok)
	WebRates(x, y, t, c, ic, cdot, 0, wdata)

def Precond(t, c, fc, jok, jcurPtr, gamma, P_data, vtemp1, vtemp2, vtemp3):
	wdata = ctypes.cast(P_data, PWebData).contents
	cvode_mem = wdata.cvode_memF
	rewt = nvecserial.NVector(wdata.rewt)
	cvodes.CVodeGetErrWeights(cvode_mem, rewt)
	
	uround = nvecserial.UNIT_ROUNDOFF
	
	fac = fc.wrmsnorm(rewt)
	r0 = 1000.0*abs(gamma)*uround*(NEQ+1)*fac
	if r0 == 0.0:
		r0 = 1.0
	
	for igy in range(wdata.ngy):
		jy = wdata.jyr[igy]
		if00 = jy*wdata.mxmp
		for igx in range(wdata.ngx):
			jx = wdata.jxr[igx]
			if0 = if00 + jx*wdata.mp
			ig = igx + igy*wdata.ngx
			for j in range(wdata.mp):
				jj = if0 + j
				save = c[jj]
				r = max([wdata.srur*abs(save),r0/rewt[jj]])
				c[jj] += r
				fac = -gamma/r
				fblock (t, c, jx, jy, vtemp1, wdata)
				for i in range(wdata.mp):
					wdata.P[ig][j][i] = (vtemp1[i] - wdata.fsave[if0+i])*fac
				c[jj] = save
	
	for ig in range(wdata.ngrp):
		cvodes.denaddI(wdata.P[ig], wdata.mp)
		ier = cvodes.denGETRF(wdata.P[ig], wdata.mp, wdata.mp, wdata.pivot[ig])
		if ier != 0:
			return 1
	
	jcurPtr.contents.value = 1
	return 0

def PrecondB(t, c, cB, fcB, jok, jcurPtr, gamma, P_data, vtemp1, vtemp2, vtemp3):
	wdata = ctypes.cast(P_data, PWebData).contents
	cvode_mem = cvodes.CVadjGetCVodeBmem(wdata.cvadj_mem)
	rewt = nvecserial.NVector(wdata.rewt)
	cvodes.CVodeGetErrWeights(cvode_mem, rewt)
	
	uround = nvecserial.UNIT_ROUNDOFF
	
	fac = fcB.wrmsnorm(rewt)
	r0 = 1000.0*abs(gamma)*uround*NEQ*fac
	if r0 == 0.0:
		r0 = 1.0
	
	for igy in range(wdata.ngy):
		jy = wdata.jyr[igy]
		if00 = jy*wdata.mxmp
		for igx in range(wdata.ngx):
			jx = wdata.jxr[igx]
			if0 = if00 + jx*wdata.mp
			ig = igx + igy*wdata.ngx
			for j in range(wdata.mp):
				jj = if0 + j
				save = c[jj]
				r = max([wdata.srur*abs(save),r0/rewt[jj]])
				c[jj] += r
				fac = -gamma/r
				fblock (t, c, jx, jy, vtemp1, wdata)
				for i in range(wdata.mp):
					wdata.P[ig][j][i] = (vtemp1[i] - wdata.fsave[if0+i])*fac
				c[jj] = save
	
	for ig in range(wdata.ngrp):
		cvodes.denaddI(wdata.P[ig], wdata.mp)
		ier = cvodes.denGETRF(wdata.P[ig], wdata.mp, wdata.mp, wdata.pivot[ig])
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
	
	for i in range(wdata.ns):
		temp = 1.0/(1.0 + 2.0*gamma*(wdata.cox[i] + wdata.coy[i]))
		beta[i] = gamma*wdata.cox[i]*temp
		beta2[i] = 2.0*beta[i]
		gam[i] = gamma*wdata.coy[i]*temp
		gam2[i] = 2.0*gam[i]
		cof1[i] = temp
	
	for jy in range(wdata.my):
		iyoff = wdata.mxns*jy
		for jx in range(wdata.mx):
			ic = iyoff + wdata.ns*jx
			v_prod(x, ic, cof1, 0, z, ic, wdata.ns)
	z[:] = [0]*len(z)
	
	for iter in range(1, ITMAX+1):
		if iter > 1:
			for jy in range(wdata.my):
				iyoff = wdata.mxns*jy
				for jx in range(wdata.mx):
					ic = iyoff + wdata.ns*jx
					if jx == 0:
						x_loc = 0
					else:
						if jx == wdata.mx-1:
							x_loc = 2
						else:
							xloc = 1
					if jy == 0:
						y_loc = 0
					else:
						if jy == wdata.my-1:
							y_loc = 2
						else:
							yloc = 1
					if (3*y_loc+x_loc) == 0:
						v_sum_prods(x, ic, beta2, 0, x, ic+wdata.ns, gam2, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 1: 
						v_sum_prods(x, ic, beta, 0, x, ic+wdata.ns, gam2, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 2: 
						v_prod(x, ic, gam2, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 3: 
						v_sum_prods(x, ic, beta2, 0, x, ic+wdata.ns, gam, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 4: 
						v_sum_prods(x, ic, beta, 0, x, ic+wdata.ns, gam, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 5: 
						v_prod(x, ic, gam, 0, x, ic+wdata.mxns, wdata.ns)
					elif (3*y_loc+x_loc) == 6: 
						v_prod(x, ic, beta2, 0, x, ic+wdata.ns, wdata.ns)
					elif (3*y_loc+x_loc) == 7: 
						v_prod(x, ic, beta, 0, x, ic+wdata.ns, wdata.ns)
					elif (3*y_loc+x_loc) == 8: 
						v_zero(x, ic, wdata.ns)
		
		for jy in range(wdata.my):
			iyoff = wdata.mxns*jy
			for jx in range(wdata.mx):
				ic = iyoff + wdata.ns*jx
				if jx == 0:
					x_loc = 0
				else:
					if jx == wdata.mx-1:
						x_loc = 2
					else:
						xloc = 1
				if jy == 0:
					y_loc = 0
				else:
					if jy == wdata.my-1:
						y_loc = 2
					else:
						yloc = 1
				if (3*y_loc+x_loc) == 0: 
				  pass
				elif (3*y_loc+x_loc) == 1: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-wdata.ns, wdata.ns)
				elif (3*y_loc+x_loc) == 2: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-wdata.ns, wdata.ns)
				elif (3*y_loc+x_loc) == 3: 
				  v_inc_by_prod(x, ic, gam, 0, x, ic-wdata.mxns, wdata.ns)
				elif (3*y_loc+x_loc) == 4: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-wdata.ns, wdata.ns)
				  v_inc_by_prod(x, ic, gam, 0, x, ic-wdata.mxns, wdata.ns)
				elif (3*y_loc+x_loc) == 5: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-wdata.ns, wdata.ns)
				  v_inc_by_prod(x, ic, gam, 0, x, ic-wdata.mxns, wdata.ns)
				elif (3*y_loc+x_loc) == 6: 
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-wdata.mxns, wdata.ns)
				elif (3*y_loc+x_loc) == 7: 
				  v_inc_by_prod(x, ic, beta, 0, x, ic-wdata.ns, wdata.ns)
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-wdata.mxns, wdata.ns)
				elif (3*y_loc+x_loc) == 8: 
				  v_inc_by_prod(x, ic, beta2, 0, x, ic-wdata.ns, wdata.ns)
				  v_inc_by_prod(x, ic, gam2, 0, x, ic-wdata.mxns, wdata.ns)
		
		z[:] = z.linearsum(1.0, 1.0, x)

def PSolve(tn, c, fc, r, z, gamma, delta, lr, P_data, vtemp):
	wdata = ctypes.cast(P_data, PWebData).contents
	
	z[:] = r
	
	GSIter(gamma, z, vtemp, wdata)
	
	iv = 0
	for jy in range(wdata.my):
		igy = wdata.jigy[jy]
		for jx in range(wdata.mx):
			igx = wdata.jigx[jx]
			ig = igx + igy*wdata.ngx
			cvodes.denGETRS(wdata.P[ig], wdata.mp, wdata.pivot[ig], z.ptrto(iv))
			iv += wdata.mp
	
	z[NEQ] = r[NEQ] + gamma*doubleIntgr(z, ISPEC, wdata)
	return 0

def PSolveB(t, c, cB, fcB, r, z, gamma, delta, lr, P_data, vtemp):
	wdata = ctypes.cast(P_data, PWebData).contents
	
	z[:] = r
	
	GSIter(-gamma, z, vtemp, wdata)
	
	iv = 0
	for jy in range(wdata.my):
		igy = wdata.jigy[jy]
		for jx in range(wdata.mx):
			igx = wdata.jigx[jx]
			ig = igx + igy*wdata.ngx
			cvodes.denGETRS(wdata.P[ig], wdata.mp, wdata.pivot[ig], z.ptrto(iv))
			iv += wdata.mp
	
	return 0

tnext = 1e-08
t = cvodes.realtype(T0)
ncheck = ctypes.c_int(0)
rewt = cvodes.NVector([0]*(NEQ+1))
wdata = WebData()
for i in range(NGRP):
	wdata.P[i] = cvodes.denalloc(NS,NS)
	wdata.pivot[i] = cvodes.denallocpiv(NS)
wdata.rewt = rewt.data
InitUserData(wdata)

c = cvodes.NVector([0]*(NEQ+1))
CInit(c, wdata)

print "\nCreate and allocate CVODE memory for forward run"
cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
wdata.cvode_memF = cvode_mem.obj
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(wdata))
cvodes.CVodeMalloc(cvode_mem, f, T0, c, cvodes.CV_SS, RTOL, ATOL)

cvodes.CVSpgmr(cvode_mem, cvodes.PREC_LEFT, 0)
cvodes.CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve, ctypes.pointer(wdata))

print "\nAllocate global memory"
cvadj_mem = cvodes.CVadjMalloc(cvode_mem, NSTEPS, cvodes.CV_HERMITE)
wdata.cvadj_mem = cvadj_mem;

print "\nForward integration"
cvodes.CVodeF(cvadj_mem, TOUT, c, ctypes.byref(t), cvodes.CV_NORMAL, ctypes.byref(ncheck))

print "\n   G = int_t int_x int_y c%d(t,x,y) dx dy dt = %f \n"%(ISPEC, c[NEQ])

cB = cvodes.NVector([0]*NEQ)

print "\nCreate and allocate CVODES memory for backward run"
cvodes.CVodeCreateB(cvadj_mem, cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeSetFdataB(cvadj_mem, ctypes.pointer(wdata))
cvodes.CVodeSetMaxNumStepsB(cvadj_mem, 1000)
cvodes.CVodeMallocB(cvadj_mem, fB, TOUT, cB, cvodes.CV_SS, RTOL, ATOL)

cvodes.CVSpgmrB(cvadj_mem, cvodes.PREC_LEFT, 0)
cvodes.CVSpilsSetPreconditionerB(cvadj_mem, PrecondB, PSolveB, ctypes.pointer(wdata))

print "\nBackward integration"
cvodes.CVodeB(cvadj_mem, T0, cB, ctypes.byref(t), cvodes.CV_NORMAL)

PrintOutput(cB, NS, MXNS, wdata)

for i in range(wdata.ngrp):
	cvodes.denfree(wdata.P[i])
	cvodes.denfreepiv(wdata.pivot[i])
