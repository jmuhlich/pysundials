#sundials_core.py is part of the PySUNDIALS package, and is released under the
#following terms and conditions.

#Copyright (c) 2007, James Dominy, Brett Olivier, Jan Hendrik Hofmeyr, Johann Rohwer
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. Neither the name of the <ORGANIZATION> nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.

"""Provides low level auxiliary functonality for use by other PySUNDIALS modules,
including dynamic library linking, configuration location and parsing, and type
mapping for the realtype.
"""

import ctypes
from ctypes import util
import os
import sys
import re
import subprocess

if os.name == "nt":
	libc = ctypes.CDLL(util.find_library('msvcrt'))
else:
	libc = ctypes.CDLL(util.find_library('c'))

try:
	libc.fdopen.argtypes = [ctypes.c_int, ctypes.c_int]
	libc.fdopen.restype = ctypes.c_void_p
	fdopen = libc.fdopen
except:
	libc._fdopen.argtypes = [ctypes.c_int, ctypes.c_int]
	libc._fdopen.restype = ctypes.c_void_p
	fdopen = libc._fdopen

libsigs = {
	"nvecserial": "N_VNew_Serial",
	"nvecparallel": None,
	"cvode": "CVodeCreate",
	"cvodes": "CVodeCreate",
	"ida": "IDACreate",
	"kinsol": "KINCreate"
}

def loadlib(libname):
	'''Links the specified library into the running python interpreter.\nlibname must be one of 'c', 'aux', 'nvecserial', 'nvecparallel', 'cvode', 'cvodes', 'ida', or 'kinsol'.'''
	p = subprocess.Popen(['sundials-config', '-m', libname, '-t', 's', '-l', 'c', '-s', 'libs'], stdout=subprocess.PIPE)
	libdir = p.communicate()[0].split()[0][2:]

	found = False
	for candidate in [os.path.join(libdir,fname) for fname in os.listdir(libdir) if fname.startswith("libsundials_"+libname)]:
		try:
			lib = ctypes.CDLL(candidate)
			lib.__getattr__(libsigs[libname])
			found = True
			break
		except OSError, e:
			pass

	if not found:
		raise OSError("%s\nCannot load shared library %s. Please check you config file and ensure the paths to the shared libraries are correct."%(e, libpaths[libname]))
	else:
		return lib

auxlibname = open(os.path.join(os.path.dirname(__file__),'auxlibname'), 'r').read()
sundials_core_aux = ctypes.CDLL(os.path.join(os.path.dirname(__file__),auxlibname))

realsize = sundials_core_aux.getsizeofrealtype()

try:
	import numpy
	numpy_imported = True
	numpy_ndarray = numpy.ndarray
	from_memory = ctypes.pythonapi.PyBuffer_FromReadWriteMemory
	from_memory.restype = ctypes.py_object
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
