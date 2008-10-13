#!/usr/bin/python

try:
	from pysundials import cvodes
except ImportError:
	import cvodes
import ctypes

class UserData (ctypes.Structure):
	_fields_ = [("p", cvodes.realtype*3)]
PUserData = ctypes.POINTER(UserData)

def f(t, y, ydot, f_data):
	data = ctypes.cast(f_data, PUserData)

	yd1 = ydot[0] = -data.contents.p[0]*y[0] + data.contents.p[1]*y[1]*y[2]
	yd3 = ydot[2] = data.contents.p[2]*y[1]*y[1]
	ydot[1] = -yd1 - yd3

	return 0

def Jac(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):
	data = ctypes.cast(jac_data, PUserData)
	
	J[0][0] = -data.contents.p[0]
	J[0][1] = data.contents.p[1]*y[2]
	J[0][2] = data.contents.p[1]*y[1]
	J[1][0] = data.contents.p[0]
	J[1][1] = -data.contents.p[1]*y[2]-2*data.contents.p[2]*y[1]
	J[1][2] = -data.contents.p[1]*y[1]
	J[2][1] = 2*data.contents.p[2]*y[1]
	
	return 0

def fQ(t, y, qdot, fQ_data):
	qdot[0] = y[2]

	return 0
 
def ewt(y, w, e_data):
	rtol = 1.0e-6
	atol = [1.0e-8, 1.0e-14, 1.0e-6]
	
	for i in range(3):
		yy = y[i]
		ww = rtol * abs(yy) + atol[i]
		if ww <= 0:
			return -1
		w[i] = 1.0/ww
	
	return 0

def fB(t, y, yB, yBdot, f_dataB):
	data = ctypes.cast(f_dataB, PUserData)
	l21 = yB[1]-yB[0]
	l32 = yB[2]-yB[1]
	y23 = y[1]*y[2]
	
	yBdot[0] = -data.contents.p[0]*l21
	yBdot[1] = data.contents.p[1]*y[2]*l21 - 2.0*data.contents.p[2]*y[1]*l32
	yBdot[2] = data.contents.p[1]*y[1]*l21 - 1.0
	
	return 0

def JacB(NB, JB, t, y, yB, fyB, jac_dataB, tmp1B, tmp2B, tmp3B):
	data = ctypes.cast(jac_dataB, PUserData)

	JB[0][0] = data.contents.p[0]
	JB[0][1] = -data.contents.p[0] 
	JB[1][0] = -data.contents.p[1]*y[2]
	JB[1][1] = data.contents.p[1]*y[2]+2.0*data.contents.p[2]*y[1]
	JB[1][2] = -2.0*data.contents.p[2]*y[1]
	JB[2][0] = -data.contents.p[1]*y[1]
	JB[2][1] = data.contents.p[1]*y[1]
	
	return 0

def fQB(t, y, yB, qBdot, fQ_dataB):
	data = ctypes.cast(fQ_dataB, PUserData)

	l21 = yB[1]-yB[0];
	l32 = yB[2]-yB[1];
	y23 = y[1]*y[2];

	qBdot[0] = y[0]*l21;
	qBdot[1] = - y23*l21;
	qBdot[2] = y[1]*y[1]*l32;
	
	return 0

def PrintOutput(tfinal, yB, qB):
	print "--------------------------------------------------------"
	print "tB0:        %12.4le"%(tfinal)
	print "dG/dp:      %12.4le %12.4le %12.4le"%(-qB[0], -qB[1], -qB[2])
	print "lambda(t0): %12.4le %12.4le %12.4le"%(yB[0], yB[1], yB[2])
	print "--------------------------------------------------------"
	print

print "\nAdjoint Sensitivity Example for Chemical Kinetics"
print "-------------------------------------------------\n"
print "ODE: dy1/dt = -p1*y1 + p2*y2*y3"
print "     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2"
print "     dy3/dt =  p3*(y2)^2\n"
print "Find dG/dp for"
print "     G = int_t0^tB0 g(t,p,y) dt"
print "     g(t,p,y) = y3\n\n"

data = UserData()
data.p[0] = 0.04
data.p[1] = 1.0e4
data.p[2] = 3.0e7

y = cvodes.NVector([1.0, 0.0, 0.0])

q = cvodes.NVector([0.0])

reltolQ = cvodes.realtype(1.0e-6)
abstolQ = cvodes.realtype(1.0e-6)

print "Create and allocate CVODES memory for forward runs"

cvode_mem = cvodes.CVodeCreate(cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeMalloc(cvode_mem, f, 0.0, y, cvodes.CV_WF, 0.0, None)

cvodes.CVodeSetEwtFn(cvode_mem, ewt, None)
cvodes.CVodeSetFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVDense(cvode_mem, 3)
cvodes.CVDenseSetJacFn(cvode_mem, Jac, ctypes.pointer(data))

cvodes.CVodeQuadMalloc(cvode_mem, fQ, q)
cvodes.CVodeSetQuadFdata(cvode_mem, ctypes.pointer(data))
cvodes.CVodeSetQuadErrCon(cvode_mem, True, cvodes.CV_SS, reltolQ, abstolQ)

time = cvodes.realtype(0.0)
ncheck = ctypes.c_int(0)
steps = 150
cvadj_mem = cvodes.CVadjMalloc(cvode_mem, steps, cvodes.CV_HERMITE)

print "Forward integration ... "
cvodes.CVodeF(cvadj_mem, 4.0e7, y, ctypes.byref(time), cvodes.CV_NORMAL, ctypes.byref(ncheck))
print "done ( nst = %ld )"%(cvodes.CVodeGetNumSteps(cvode_mem))
cvodes.CVodeGetQuad(cvode_mem, 4.0e7, q)
print "--------------------------------------------------------"
print "G:          %12.4le "%(q[0])
print "--------------------------------------------------------\n"

#  /* Test check point linked list 
#     (uncomment next block to print check point information) */
#  
#  /*
#  printf("\nList of Check Points (ncheck = %d)\n\n", ncheck);
#  ckpnt = (CVadjCheckPointRec *) malloc ( (ncheck+1)*sizeof(CVadjCheckPointRec));
#  CVadjGetCheckPointsInfo(cvadj_mem, ckpnt);
#  for (i=0;i<=ncheck;i++) {
#    printf("Address:       %p\n",ckpnt[i].my_addr);
#    printf("Next:          %p\n",ckpnt[i].next_addr);
#    printf("Time interval: %le  %le\n",ckpnt[i].t0, ckpnt[i].t1);
#    printf("Step number:   %ld\n",ckpnt[i].nstep);
#    printf("Order:         %d\n",ckpnt[i].order);
#    printf("Step size:     %le\n",ckpnt[i].step);
#    printf("\n");
#  }
#  */

yB = cvodes.NVector([0.0, 0.0, 0.0])
qB = cvodes.NVector([0.0, 0.0, 0.0])
reltolB = 1.0e-6
abstolB = cvodes.realtype(1.0e-8)
abstolQB = cvodes.realtype(1.0e-6)

print "Create and allocate CVODES memory for backward run"

cvodes.CVodeCreateB(cvadj_mem, cvodes.CV_BDF, cvodes.CV_NEWTON)
cvodes.CVodeMallocB(cvadj_mem, fB, 4.0e7, yB, cvodes.CV_SS, reltolB, abstolB)
cvodes.CVodeSetFdataB(cvadj_mem, ctypes.pointer(data))
cvodes.CVDenseB(cvadj_mem, 3)
cvodes.CVDenseSetJacFnB(cvadj_mem, JacB, ctypes.pointer(data))
cvodes.CVodeQuadMallocB(cvadj_mem, fQB, qB)
cvodes.CVodeSetQuadFdataB(cvadj_mem, ctypes.pointer(data))
cvodes.CVodeSetQuadErrConB(cvadj_mem, True, cvodes.CV_SS, reltolB, abstolQB)

print "Backward integration ...",
cvodes.CVodeB(cvadj_mem, 0.0, yB, ctypes.byref(time), cvodes.CV_NORMAL)
print "done ( nst = %ld )"%(cvodes.CVodeGetNumSteps(cvodes.CVadjGetCVodeBmem(cvadj_mem)))
cvodes.CVodeGetQuadB(cvadj_mem, qB)
PrintOutput(4.0e7, yB, qB)

yB = cvodes.NVector([0.0, 0.0, 0.0])
qB = cvodes.NVector([0.0, 0.0, 0.0])

print "Re-initialize CVODES memory for backward run"

cvodes.CVodeReInitB(cvadj_mem, fB, 50.0, yB, cvodes.CV_SS, reltolB, abstolB)
cvodes.CVodeQuadReInitB(cvadj_mem, fQB, qB) 

print "Backward integration ...",
cvodes.CVodeB(cvadj_mem, 0.0, yB, ctypes.byref(time), cvodes.CV_NORMAL)
print "done ( nst = %ld )"%(cvodes.CVodeGetNumSteps(cvodes.CVadjGetCVodeBmem(cvadj_mem)))
cvodes.CVodeGetQuadB(cvadj_mem, qB)
PrintOutput(50.0, yB, qB)
