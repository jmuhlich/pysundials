#nvecserial.py is part of the PySUNDIALS package, and is released under the
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

"""Python bindings for the serial NVector type

This module exposes two classes NVector and NVectorArray"""

import ctypes
import sundials_core
import math
import sys

numpy_imported = sundials_core.numpy_imported
realtype = sundials_core.realtype
UNIT_ROUNDOFF = sundials_core.UNIT_ROUNDOFF

if numpy_imported:
	numpy_ndarray = sundials_core.numpy_ndarray
	numpyrealtype = sundials_core.numpyrealtype

nvecserial = sundials_core.loadlib("nvecserial")

class _N_VectorContent_Serial(ctypes.Structure):
	_fields_ = [("length", ctypes.c_long), ("owndata", ctypes.c_int), ("data", ctypes.POINTER(realtype))]

class _NVector(ctypes.Structure):
	_fields_ = [("content", ctypes.POINTER(_N_VectorContent_Serial)), ("ops", ctypes.c_void_p)]
PVector = ctypes.POINTER(_NVector)

nvecserial.N_VNew_Serial.restype = PVector

nvecserial.N_VLinearSum_Serial.argtypes = [realtype, PVector, realtype, PVector, PVector]
nvecserial.N_VLinearSum_Serial.restype = None

nvecserial.N_VProd_Serial.argtypes = [PVector, PVector, PVector]
nvecserial.N_VProd_Serial.restype = None

nvecserial.N_VDiv_Serial.argtypes = [PVector, PVector, PVector]
nvecserial.N_VDiv_Serial.restype = None

nvecserial.N_VScale_Serial.argtypes = [realtype, PVector, PVector]
nvecserial.N_VScale_Serial.restype = None

nvecserial.N_VAbs_Serial.argtypes = [PVector, PVector]
nvecserial.N_VAbs_Serial.restype = None

nvecserial.N_VAddConst_Serial.argtypes = [PVector, realtype, PVector]
nvecserial.N_VAddConst_Serial.restype = None

nvecserial.N_VDotProd_Serial.argtypes = [PVector, PVector]
nvecserial.N_VDotProd_Serial.restype = realtype

nvecserial.N_VWrmsNorm_Serial.argtypes = [PVector, PVector]
nvecserial.N_VWrmsNorm_Serial.restype = realtype

nvecserial.N_VWrmsNormMask_Serial.argtypes = [PVector, PVector, PVector]
nvecserial.N_VWrmsNormMask_Serial.restype = realtype

nvecserial.N_VWL2Norm_Serial.argtypes = [PVector, PVector]
nvecserial.N_VWL2Norm_Serial.restype = realtype

nvecserial.N_VL1Norm_Serial.argtypes = [PVector]
nvecserial.N_VL1Norm_Serial.restype = realtype

nvecserial.N_VCompare_Serial.argtypes = [realtype, PVector, PVector]
nvecserial.N_VCompare_Serial.restype = None

nvecserial.N_VInvTest_Serial.argtypes = [PVector, PVector]
nvecserial.N_VInvTest_Serial.restype = ctypes.c_int

nvecserial.N_VConstrMask_Serial.argtypes = [PVector, PVector, PVector]
nvecserial.N_VConstrMask_Serial.restype = ctypes.c_int

nvecserial.N_VMinQuotient_Serial.argtypes = [PVector, PVector]
nvecserial.N_VMinQuotient_Serial.restype = realtype

class NVector(object):
	"""The NVector object provides a convenient wrapper around the SUNDIALS NVector structure. It can be indexed or sliced like any other python sequenece, and operations on NVectors return new object, not modifying the original operands. A complete set of vector operations are implemented, as operator overloads where intuitive."""
	def __init__(self, vector):
		"""NVector.__init__(self, vector) -> NVector; vector may be an NVector, numpy ndarray of appropriate type, or any python seqenece of real numbers\n\nAdditionally, vector may be of the special type LP__NVector, which is a ctypes pointer to a C _NVector structure, in which case new memory is not allocated for the NVector, but only the wrapping python NVector object is created around the same underlying _NVector. This second use exists primarily for internal PySUNDIALS use. Use with care."""
		if type(vector) in [list, NVector]:
			self.length = len(vector)
			self.data = nvecserial.N_VNew_Serial(len(vector))
			self.copy = False
			for v in range(len(vector)):
				self.data.contents.content.contents.data[v] = vector[v]
		elif type(vector).__name__ == "LP__NVector":
			self.length = vector.contents.content.contents.length
			self.data = vector
			self.copy = True
		elif numpy_imported and type(vector) == numpy_ndarray:
			if vector.ndim > 1:
				raise TypeError("Cannot create NVector from ndarray of dimension greater than 1")
			if vector.dtype != numpyrealtype:
				raise TypeError("Cannot create NVector from ndarray of non-matching dtype (%s)"%vector.dtype)
			self.length = len(vector)
			self.data = nvecserial.N_VNew_Serial(len(vector))
			self.copy = False
			ctypes.memmove(self.addressof(), vector.ctypes.data, ctypes.sizeof(realtype)* self.length)
		else:
			raise TypeError("Cannot create NVector from type %s"%(type(vector).__name__))
	
	def __getitem__(self, index):
		"""x.__getitem__(y) <==> x[y]"""
		if (index < 0) or (index >= self.length):
			raise IndexError("Vector index out of bounds")
		return self.data.contents.content.contents.data[index]

	def __setitem__(self, index, value):
		"""x.__setitem__(i, y) <==> x[i] = y"""
		if (index < 0) or (index >= self.length):
			raise IndexError("Vector index out of bounds")
		self.data.contents.content.contents.data[index] = value

	def __getslice__(self, i, j):
		"""x.__getslice__(i, j) <==> x[i:j]"""
		if (i < 0):
			i += self.length
		if (j < 0):
			j += self.length
		ret = NVector([0]*(j-i))
		x = 0
		while x+i < j:
			ret[x] = self[x+i]
			x += 1
		return ret

	def __setslice__(self, i, j, y):
		"""x.__setslice__(i, j, y) <==> x[i:j] = y"""
		if (i < 0):
			i += self.length
		if (j < 0):
			j += self.length
		if (j == sys.maxint):
			j = self.length
		x = i
		while x < j:
			self[x] = y[x-i]
			x += 1

	def __len__(self):
		"""x.__len__() <==> len(x)"""
		return self.length

	def __del__(self):
		"""NVector destructor: Frees any memory allocated in C"""
		if not self.copy:
			nvecserial.N_VDestroy_Serial(self.data)
		
	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		s = "["
		for e in range(self.length-1):
			s += repr(self.data.contents.content.contents.data[e]) + ", "
		s += repr(self.data.contents.content.contents.data[self.length-1]) + "]"
		return s
	
	def __neg__(self): #unary element wise negation
		"""x.__neg__() <==> -x (Unary element wise negation)"""
		ret = NVector(self)
		nvecserial.N_VScale_Serial(-1, self.data, ret.data)
		return ret

	def __abs__(self): #element wise absolute value
		"""x.__abs__() <==> abs(x) (Element wise absolute value)"""
		ret = NVector(self)
		nvecserial.N_VAbs_Serial(self.data, ret.data)
		return ret
	
	def __add__(self, v): #scalar and vector addition
		"""x.__add__(y) <==> x+y (y is scalar: element wise addition of y)\nx.__add__(y) <==> x+y (y is a sequence: vector addition of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VAddConst_Serial(self.data, v, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot add two vectors of different lengths")
			for i in range(ret.length):
				ret[i] = self[i] + v[i]
		else:
			raise TypeError("Cannot add NVector and %s"%(type(v).__name__))
		return ret
	
	def __radd__(self, v): #scalar and vector addition
		"""x.__add__(y) <==> x+y (y is scalar: element wise addition of y)\nx.__add__(y) <==> x+y (y is a sequence: vector addition of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VAddConst_Serial(self.data, v, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot add two vectors of different lengths")
			for i in range(ret.length):
				ret[i] = self[i] + v[i]
		else:
			raise TypeError("Cannot add NVector and %s"%(type(v).__name__))
		return ret
	
	def __sub__(self, v): #scalar and vector subtraction
		"""x.__sub__(y) <==> x-y (y is scalar: element wise subtraction of y)\nx.__sub__(y) <==> x-y (y is a sequence: vector subtraction of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VAddConst_Serial(self.data, -v, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot subtract two vectors of different lengths")
			for i in range(ret.length):
				ret[i] = self[i] - v[i]
		else:
			raise TypeError("Cannot subtract %s from NVector"%(type(v).__name__))
		return ret

	def __rsub__(self, v): #vector subtraction
		"""x.__sub__(y) <==> x-y (y is scalar: element wise subtraction of y)\nx.__sub__(y) <==> x-y (y is a sequence: vector subtraction of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot subtract two vectors of different lengths")
			for i in range(ret.length):
				ret[i] = v[i] - self[i]
		else:
			raise TypeError("Cannot subtract NVector from %s"%(type(v).__name__))
		return ret

	def __mul__(self, v): #scalar and vector multiplication
		"""x.__mul__(y) <==> x*y (y is scalar: element wise multiplication of y [Scaling])\nx.__mul__(y) <==> x*y (y is a sequence: vector multiplication of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VScale_Serial(v, self.data, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot multiply two vectors of different lengths")
			if type(v).__name__ in ['tuple', 'list']:
				vec = NVector(v)
			else:
				vec = v
			nvecserial.N_VProd_Serial(self.data, vec.data, ret.data)
		else:
			raise TypeError("Cannot multiply NVector and %s"%(type(v).__name__))
		return ret

	def __rmul__(self, v): #scalar and vector multiplication
		"""x.__mul__(y) <==> x*y (y is scalar: element wise multiplication of y [Scaling])\nx.__mul__(y) <==> x*y (y is a sequence: vector multiplication of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VScale_Serial(v, self.data, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot multiply two vectors of different lengths")
			if type(v).__name__ in ['tuple', 'list']:
				vec = NVector(v)
			else:
				vec = v
			nvecserial.N_VProd_Serial(self.data, vec.data, ret.data)
		else:
			raise TypeError("Cannot multiply NVector and %s"%(type(v).__name__))
		return ret

	def __div__(self, v): #scalar and vector division
		"""x.__div__(y) <==> x/y (y is scalar: element wise division of y [Scaling])\nx.__div__(y) <==> x/y (y is a sequence: vector division of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['float', 'int']:
			nvecserial.N_VScale_Serial(1.0/v, self.data, ret.data)
		elif type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot divide two vectors of different lengths")
			if type(v).__name__ in ['tuple', 'list']:
				vec = NVector(v)
			else:
				vec = v
			nvecserial.N_VDiv_Serial(self.data, vec.data, ret.data)
		else:
			raise TypeError("Cannot divide NVector by %s"%(type(v).__name__))
		return ret

	def __rdiv__(self, v): #vector division
		"""x.__div__(y) <==> x/y (y is scalar: element wise division of y [Scaling])\nx.__div__(y) <==> x/y (y is a sequence: vector division of y)"""
		ret = NVector(self)
		if type(v).__name__ in ['NVector', 'tuple', 'list']:
			if len(v) != self.length:
				raise TypeError("Cannot divide two vectors of different lengths")
			if type(v).__name__ in ['tuple', 'list']:
				vec = NVector(v)
			else:
				vec = v
			nvecserial.N_VDiv_Serial(vec.data, self.data, ret.data)
		else:
			raise TypeError("Cannot divide %s by NVector"%(type(v).__name__))
		return ret

	def addressof(self, index = 0):
		"""Returns the address of a particular realtype of index 'index' from within the NVector's actual data array. Useful for passing a pointer to a partiular portion of the NVector."""
		return ctypes.addressof(self.data.contents.content.contents.data.contents)+(index * ctypes.sizeof(realtype))
	
	def ptrto(self, index = 0):
		return ctypes.pointer(realtype.from_address(self.addressof(index)))

	def inverted(self): #element wise 1/x
		"""x.__invert__() <==> 1/x (Element wise 1/x)"""
		ret = NVector(self)
		if nvecserial.N_VInvTest_Serial(self.data, ret.data) == 1:
			return ret
		else:
			raise ZeroDivisionError("Vectors containing elements of value zero cannot be inverted")

	def linearsum(self, a, b, y):
		"""x.linearsum(self, a, b, y) <==> a*x + b*y"""
		if type(y).__name__ in ['tuple', 'list']:
			v = NVector(y)
		elif type(y).__name__ == 'NVector':
			v = y
		else:
			raise TypeError('y must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		ret = NVector(self)
		nvecserial.N_VLinearSum_Serial(a, self.data, b, v.data, ret.data)
		return ret

	def dotproduct(self, y):
		"""x.dotproduct(self, y) <==> the dot product of x and y"""
		if type(y).__name__ in ['tuple', 'list']:
			v = NVector(y)
		elif type(y).__name__ == 'NVector':
			v = y
		else:
			raise TypeError('y must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		if len(v) != self.length:
			raise TypeError("Cannot determine dot product of two vectors of different lengths")
		return nvecserial.N_VDotProd_Serial(self.data, v.data)
	
	def wrmsnorm(self, w):
		"""x.wrmsnorm(self, w) <==> root mean square of x weighted by v"""
		if type(w).__name__ in ['tuple', 'list']:
			v = NVector(w)
		elif type(w).__name__ == 'NVector':
			v = w
		else:
			raise TypeError('w must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		return nvecserial.N_VWrmsNorm_Serial(self.data, v.data)

	def wrmsnormmask(self, w, mask): 
		"""x.wrmsnormmask(self, w, mask) <==> root mean sqiare of x weighted by w and masked by mask"""
		if type(w).__name__ in ['tuple', 'list']:
			v = NVector(w)
		elif type(w).__name__ == 'NVector':
			v = w
		else:
			raise TypeError('w must be a sequence')
		if type(mask).__name__ in ['tuple', 'list']:
			m = NVector(mask)
		elif type(mask).__name__ == 'NVector':
			m = mask
		else:
			raise TypeError('mask must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		if self.length != len(m):
			raise TypeError("Vectors must be the same length")
		return nvecserial.N_VWrmsNormMask_Serial(self.data, v.data, m.data)

	def wl2norm(self, w = None):
		"""x.wl2norm(self, w = None) <==> euclidean L2 norm of x, weighted by vector w; If w is not given, weighting is equal"""
		if w is None:
			w = NVector([1]*self.length)
		if type(w).__name__ in ['tuple', 'list']:
			wv = NVector(w)
		else:
			wv = w
		if self.length != len(wv):
			raise TypeError("Arguments must be sequences of the same length")
		return nvecserial.N_VWL2Norm_Serial(self.data, wv.data)
	
	def l1norm(self): 
		"""x.l1norm(self) <==> sum of absolute values of the elements of x"""
		return nvecserial.N_VL1Norm_Serial(self.data)
	
	def compare(self, v): 
		"""x.compare(self, v) <==> m; where m is a vector such that x[i] => v => m[i] = 1, else m[i] = 0"""
		if type(v).__name__ not in ['float', 'int']:
			raise TypeError("Cannot compare elements of an NVector to %s"%(type(v).__name__))
		ret = NVector(self)
		nvecserial.N_VCompare_Serial(v, self.data, ret.data)
		return ret

	def constrain(self, c):
		"""x.constrain(self, c) <==> True if x all elements of x pass their respective contraint tests, otherwise m; where m is a vector such that
    m[i] = 1.0 if constraint test fails for x[i]
    m[i] = 0.0 if constraint test passes for x[i]
 where the constraint tests are as follows:
    If c[i] = +2.0, then x[i] must be >  0.0.
    If c[i] = +1.0, then x[i] must be >= 0.0.
    If c[i] = -1.0, then x[i] must be <= 0.0.
    If c[i] = -2.0, then x[i] must be <  0.0."""
		if type(c).__name__ in ['tuple', 'list']:
			v = NVector(c)
		elif type(c).__name__ == 'NVector':
			v = c
		else:
			raise TypeError('c must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		ret = NVector(self)
		if nvecserial.N_VConstrMask_Serial(v.data, self.data, ret.data) == 1:
			return True
		else:
			return ret

	def minquotient(self, d):
		"""x.minquotient(self, d) <=> min(x/d); ignoring 0 value elements in d. If all elements of d are 0, BIG_REAL is returned."""
		if type(d).__name__ in ['tuple', 'list']:
			v = NVector(d)
		elif type(c).__name__ == 'NVector':
			v = d
		else:
			raise TypeError('d must be a sequence')
		if self.length != len(v):
			raise TypeError("Vectors must be the same length")
		return nvecserial.N_VMinQuotient_Serial(self.data, v.data)

class NVectorArray(object):
	"""A wrapper around the 'array of NVectors' structure used for sensitivity analysis. Individual elements of the array are returned as python NVector objects."""
	def __init__(self, vector_init):
		"""NVectorArray.__init__(self, vector_init) -> NVectorArray; vector_init must be a python seqence of python sequences of real numbers.\n\nNB: if given a sequence of NVector objects, the internal elements of the array returned are copies of the original NVector objects."""
		self.length = len(vector_init)
		self.data = []
		for v in vector_init:
			self.data.append(NVector(v))
		self.cdata = (PVector*self.length)()
		for i in range(self.length):
			self.cdata[i] = self.data[i].data
	
	def __len__(self):
		"""x.__len__() <==> len(x)"""
		return self.length
	
	def __getitem__(self, index):
		"""x.__getitem__(y) <==> x[y]; x[y] is an NVector object"""
		return self.data[index]
	
	def __repr__(self):
		"""x.__repr__() <==> repr(x)"""
		ret = "["
		for row in range(self.length):
			ret += self.data[row].__repr__()
			if row < self.length-1:
				ret += "\n "
			else:
				ret += "]"
		return ret
