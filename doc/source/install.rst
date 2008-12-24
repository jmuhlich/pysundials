Installation
------------

Prerequisites
+++++++++++++

PySUNDIALS requires the following

 * a working copy of SUNDIALS installed with shared libraries compiled
 * Python v2.4 with the ``ctypes`` module, OR Python v2.5 or higher
 
The python packages ``numpy`` and ``scipy`` are recommended but not necessary.

From Source
+++++++++++

The latest release version of PySUNDIALS can be download at
http://sourceforge.net/projects/pysundials, or alternatively, the very latest
(potentially non-functional) version via anonymous svn using the command::

        $ svn co https://pysundials.svn.sourceforge.net/svnroot/pysundials pysundials

On Linux/BSD or other POSIX
...........................

        * Download and untar the complete SUNDIALS suite.
        * ``$ ./configure --enable-shared``
        * ``$ make && make install``
        * Download and untar PySUNDIALS
        * Change to the directory where you unpacked PySUNDIALS
        * ``$ python setup.py install``

On Windows using MSys/MinGW
...........................

        * Download and untar the complete SUNDIALS suite somewhere inside MSys.
        * Do *not* run aclocal, autoconf, or autoheader, or the makefiles will break!
        * ``$ ./configure --enable-shared``
        * ``$ make && make install``
        * Download and untar PySUNDIALS
        * Change to the directory where you unpacked PySUNDIALS
        * ``$ python setup.py -c mingw32 install``

If you receive a compile time error when executing the final command, which
complains about missing ``sundials_conf.h``, or other missing ``.h`` files, set
the ``CPATH`` environment variable to point to the directory containing the
sundials include directories, for example::

        $ CPATH=/local/include python setup.py -c mingw32 install

Note that using MSys presents unique problems being a hybrid environment. We
recommend using MSys only to compile and install (Py)SUNDIALS, not for use as
an environment in which to run PySUNDIALS. Having followed the above
instructions, the SUNDIALS shared libs will be in ``%MSYSROOT%/usr/local/lib/``,
and end with ``-<digit>.exe``.  On older versions of MSys/MinGW, the ``.exe``
extension may be left off.

Binary Installation on Windows
++++++++++++++++++++++++++++++

Configuration
+++++++++++++

PySUNDIALS uses ``ctypes`` to link the required SUNDIALS libraries directly into
the running python process. In order to do this, it needs to know where to find
those shared libraries. On Linux/BSD systems, this is usually auto detected, as
library locations are standard, however on windows systems, the final location
and even naming convention of the shared library files is compiler dependent.
If PySUNDIALS cannot find the SUNDIALS libraries, please locate them yourself,
and specify their locations in a file named ``~/.pysundials/config`` (Linux/BSD),
or ``%HOMEPATH%\pysundials\config`` (Windows). If PySUNDIALS cannot find this
configuration file in your home directory, it will seek it in the same
directory it was installed in, i.e. ``$PYTHONROOT/site-packages/pysundials/config``.

In the 'config' file, anything following a hash (including the hash itself) is
considered a comment. Each line specifies the location of a required library in
the form::

        library = path

Where *library* is one of ``c``, ``aux``, ``nvecserial``, ``nvecparallel``, ``cvode``,
``cvodes``, ``ida``, or ``kinsol``, and *path* is the complete path of the
appropriate library file. ``c``, ``aux``, and ``nvecparallel`` are optional; the
first two being generally autodetected, and the second only required if you
will be using PySUNDIALS in parallel with MPI.

You can find a sample config file for both posix and mingw32 systems in the doc
subdirectory of the PySUNDIALS distribution.
