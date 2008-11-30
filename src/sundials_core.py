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
	'''Links the specified library into the running python interpreter.\nlibname must be one of 'c', 'aux', 'nvecserial', 'nvecparallel', 'cvode', 'cvodes', 'ida', or 'kinsol'.'''
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
