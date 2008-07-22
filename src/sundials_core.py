import ctypes
from ctypes import util
import os
import sys
import re

if os.name == "nt":
	clib = util.find_library('msvcrt')
else:
	clib = util.find_library('c')

auxlibname = open(os.path.dirname(__file__)+'/auxlibname', 'r').read()
libpaths = {
	"c": clib,
	"aux": os.path.dirname(__file__)+'/'+auxlibname,
	"nvecserial": None,
	"nvecparallel": None,
	"cvode": None,
	"cvodes": None,
	"ida": None,
	"kinsol": None
}

def loadlib(libname):
	try:
		lib = ctypes.CDLL(libpaths[libname])
	except OSError, e:
		raise OSError("%s\nCannot load shared library %s. Please check you config file and ensure the paths to the shared libraries are correct."%(e, libpaths[libname]))
	return lib

try:
	if os.sys.platform == "win32":
		homedir = os.getenv("HOMEDRIVE")+'/'+os.getenv("HOMEPATH")
		f = open(homedir+"/pysundials/config", "r")
	else:
		homedir = os.getenv("HOME")
		f = open(homedir+"/.pysundials/config", "r")
except IOError:
	try:
		f = open(os.path.dirname(__file__)+"/config", "r")
	except IOError:
		raise IOError("Couldn't find pysundials config file at %s. Please create one, and point it to the appropriate library files."%(os.path.dirname(__file__)+"/config"))
	
item = re.compile('(^.*)(#.*)?$')
for line in f:
	m = item.match(line)
	if m.group(1) != '':
		lib, path = tuple([x.strip() for x in m.group(1).split('=')])
		libpaths[lib] = path
f.close()

libc = loadlib("c")
try:
	libc.fdopen.argtypes = [ctypes.c_int, ctypes.c_int]
	libc.fdopen.restype = ctypes.c_void_p
	fdopen = libc.fdopen
except:
	libc._fdopen.argtypes = [ctypes.c_int, ctypes.c_int]
	libc._fdopen.restype = ctypes.c_void_p
	fdopen = libc._fdopen

sundials_core_aux = loadlib("aux")

realsize = sundials_core_aux.getsizeofrealtype()
#copyarray = sundials_core_aux.copyarray

try:
	import numpy
	numpy_imported = True
	numpy_ndarray = numpy.ndarray
except:
	numpy_imported = False

if realsize > ctypes.sizeof(ctypes.c_double):
	raise AssertionError("SUNDIALS ERROR: size of realtype is extended, which cannot be represented in python!")
elif realsize == ctypes.sizeof(ctypes.c_double):
	realtype = ctypes.c_double
	UNIT_ROUNDOFF = 1e-9
	if numpy_imported:
		numpyrealtype = 'float64'
else:
	realtype = ctypes.c_float
	UNIT_ROUNDOFF = 1e-6
	if numpy_imported:
		numpyrealtype = 'float32'

del realsize
