import os
import sys
from distutils.core import setup
from distutils import ccompiler

def usage(errmesg):
	sys.stderr.write("""%s
Usage: %s [--compiler=<compiler>] [-I<include dir>] [<command>]
	The short form of --compiler (-c <compiler>) is also valid
	<command> is one of the usual distutils setup commands
"""%(errmesg, sys.argv[0]))
	sys.exit(1)
compiler_type = None
include_dirs = []

i = 1
while i < len(sys.argv)-1:
	arg = sys.argv[i]
	if arg.startswith("--compiler="):
		compiler_type = arg[11:]
		if compiler_type not in ccompiler.compiler_class.keys():
			usage('%s is not a valid compiler, choose from %r'%(compiler_type, ccompiler.compiler_class.keys()))
		del sys.argv[i]
	elif arg == '-c':
		del sys.argv[i]
		compiler_type = sys.argv[i]
		if compiler_type not in ccompiler.compiler_class.keys():
			usage('%s is not a valid compiler, choose from %r'%(compiler_type, ccompiler.compiler_class.keys()))
		del sys.argv[i]
	elif arg.startswith("-I"):
		include_dirs.append(arg[2:])
		del sys.argv[i]
	else:
		i += 1

compiler = ccompiler.new_compiler(None, compiler_type)
for include in include_dirs:
	compiler.add_include_dir(include)

if compiler_type == "mingw32":
	compiler.compile(['src/sundials_core_aux.c'])
	compiler.link_shared_lib(['src/sundials_core_aux.o'], 'libsundials_core_aux', 'build/lib/pysundials')
	libfname = compiler.library_filename('libsundials_core_aux', 'shared')
else:
	compiler.compile(['src/sundials_core_aux.c'], extra_postargs=['-fPIC'])
	compiler.link_shared_lib(['src/sundials_core_aux.o'], 'sundials_core_aux', 'build/lib/pysundials', None, None, None, None, 0, None, ['-fPIC'])
	libfname = compiler.library_filename('sundials_core_aux', 'shared')

f = open('src/auxlibname', "w")
f.write(libfname)
f.close()

setup(
	 name = "pysundials",
	 version = "2.3.0-b6",
	 description = "Python wrappers for the SUite of Non-linear DIfferential and ALgebraic solvers (SUNDIALS).",
	 maintainer = 'James Dominy',
	 maintainer_email = 'james@sun.ac.za',
	 url = 'http://www.sourceforge.net/projects/pysundials',
	 package_dir = {'pysundials': 'src'},
	 packages = ['pysundials'],
	 package_data = {'pysundials': ['auxlibname']},
	 classifiers = [
		'Development Status :: 4 - Beta', 
		'Environment :: Console', 
		'Intended Audience :: Science/Research', 
		'License :: OSI Approved :: MIT License',
		'Operating System :: OS Independent',
		'Programming Language :: Python',
		'Programming Language :: C',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Scientific/Engineering :: Mathematics',
		'Topic :: Software Development :: Libraries :: Python Modules'
	]
)

os.remove('src/auxlibname')
