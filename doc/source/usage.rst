Using PySUNDIALS
----------------

Importing PySUNDIALS modules
++++++++++++++++++++++++++++

There are five PySUNDIALS modules available for general use; One for each of the
SUNDIALS modules, namely ``cvode``, ``cvodes``, ``ida``, and ``kinsol``. The
fifth is the ``nvecserial`` module, which you may wish to use in conjuction with
CVODES for it's convenience functions. To import one of these modules use::

        from pysundials import module_name

where *module_name* is the name of the module you wish to use.

NVectors
++++++++

The fundamental data type of SUNDIALS, and hence PySUNDIALS is the NVector.
PySUNDIALS implements the NVector class in a manner that closely resembles a
python list. Creating an NVector is a simple case of class instantiation::

        >>> from pysundials import nvecserial
        >>> v = nvecserial.NVector([1, 2, 3])
        >>> v
        [1.0, 2.0, 3.0]

When instantiating an NVector, simply pass a sequence (tuple, list, numpy array,
another NVector, etc...) to the constuctor. The contents of the sequence
determine the length of the NVector (which remains immutable for its lifetime)
and the initial value of the NVector. NVector objects can be subscripted or
sliced like normal python lists.::

        >>> v[1]
        2.0
        >>> v[2] = 4
        >>> v
        [1.0, 2.0, 4.0]
        >>> v[0:2]
        [1.0, 2.0]
        >>> v[0:2] = (-1, -2)
        >>> v
        [-1.0, -2.0, 4.0]

Operators act intuitively on NVectors too,

        * ``+`` and ``-`` perform scalar or vector addition or subtraction
          respectively, depending on operands.::

                >>> w = nvecserial.NVector([1,2,-4])
                >>> v+1
                [0.0, -1.0, 5.0]
                >>> v+w
                [0.0, 0.0, 0.0]
                >>> v-w
                [-2.0, -4.0, 8.0]

        * ``\*`` performs scalar or element-wise multiplication depending on
          operands.::

                >>> v*2
                [-2.0, -4.0, 8.0]
                >>> v*w
                [-1.0, -4.0, -16.0]

        * ``/`` performs scalar or element-wise division depending on opreands.::

                >>> v/2
                [-0.5, -1.0, 2.0]
                >>> v/w
                [-1.0, -1.0, -1.0]
                >>> v/v
                [1.0, 1.0, 1.0]

        * and a variety of object methods perform more complex operations
          including dot product and various norms. See the reference section for a
          complete list

An NVector object can be used as a numpy array by using its ``asarray`` method.
Note how changes to the array affect the NVector and *vice versa*.::

        >>> import numpy
        >>> a = v.asarray()
        >>> a
        array([-1., -2.,  4.])
        >>> a[0] = 0
        >>> a
        array([ 0., -2.,  4.])
        >>> v
        [0.0, -2.0, 4.0]
        >>> v[1] = 0
        >>> a
        array([ 0.,  0.,  4.])

The NVector class is exported to each of the main PySUNDIALS modules, so there
is rarely a need to import ``nvecserial``.::

        >>> from pysundials import cvode
        >>> v = cvode.NVector([1,2,3])
        >>> v
        [1.0, 2.0, 3.0]

CVODE
+++++

Programs using CVODE will generally conform to a certain skeleton layout.

        #. Import the covde module::

                from pysundials import cvode
                       
        #. Define your right-hand side function. This function must take exactly
           four parameter. The first parameter will be the current value of the
           indepedent variable (usually time). The second paramter will be an
           Nvector containing the current values of the dependent variables. The
           third parameter is an NVector whose elements must be filled with the new
           values of the dependent variables. The fourth parameter is a pointer to
           any arbitrary user data you may have specified, otherwise None. This
           function essentially defines your ODE system::

                def f(t, y, ydot, f_data):
                        yd1 = ydot[0] = -0.04*y[0] + 1.0e4*y[1]*y[2]
                        yd3 = ydot[2] = 3.0e7*y[1]*y[1]
                        ydot[1] = -yd1 - yd3
                        return 0

        #. Define any optional functions such as Jacobian, error weight and/or root
           finding functions.::

                def g(t, y, gout, g_data):
                        gout[0] = y[0] - 0.0001
                        gout[1] = y[2] - 0.01
                        return 0

                def Jac(N, J, t, y, fy, jac_data, tmp1, tmp2, tmp3):
                        J[0][0] = -0.04
                        J[0][1] = 1.0e4*y[2]
                        J[0][2] = 1.0e4*y[1]
                        J[1][0] = 0.04
                        J[1][1] = -1.0e4*y[2]-6.0e7*y[1]
                        J[1][2] = -1.0e4*y[1]
                        J[2][1] = 6.0e7*y[1]
                        return 0

        #. Initialise an NVector with the initial conditions.::

                y = cvode.NVector([1.0, 0.0, 0.0])

        #. Create a CVODE object.::

                cvode_mem = cvode.CVodeCreate(cvode.CV_BDF, cvode.CV_NEWTON)
                        
        #. Allocate integrator memory, set the initial value of the independent
           variable, and set tolerances.::

                abstol = cvode.NVector([1.0e-8, 1.0e-14, 1.0e-6])
                reltol = cvode.realtype(1.0e-4)
                cvode.CVodeMalloc(cvode_mem, f, 0.0, y, cvode.CV_SV, reltol, abstol)

        #. Set any optional inputs using ``CVSet*()``.
        #. Choose a linear solver and set the problem size, i.e. number of
        variables.::

                cvode.CVDense(cvode_mem, 3)

        #. Set any optional linear solver inputs.::

                cvode.CVDenseSetJacFn(cvode_mem, Jac, None)

        #. Optionally initialise root finding.::

                cvode.CVodeRootInit(cvode_mem, 2, g, None)

        #. Advance the solution in time, calling ``cvode.CVode`` for each
           desired output time step. Each call to ``cvode.CVode`` specifies the
           desired time for the next stop (``tout``) and the current conditions
           (``y``). On return, ``y`` will contain the new conditions, and ``t``
           will contain the time at which the integrator stopped. ``t`` may be
           different from ``tout`` if roots are found, or errors encountered. The
           last parameter specifies how CVODE should step. See the SUNDIALS
           documentation for more details.::

                tout = 0.4
                while tout < 0.4*(10**12):
                        flag = cvode.CVode(cvode_mem, tout, y, ctypes.byref(t), cvode.CV_NORMAL)
                        print (t, y)
                        if flag == cvode.CV_ROOT_RETURN:
                                rootsfound = cvode.CVodeGetRootInfo(cvode_mem, 2)
                                print rootsfound
                        elseif flag == cvode.CV_SUCCESS:
                                tout *= 10
                        else:
                                break

CVODES
++++++

IDA
+++

KINSOL
++++++
