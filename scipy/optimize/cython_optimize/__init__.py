"""
Cython Optimize Zeros API
=========================
The following four functions in ``scipy.optimize`` can be accessed directly
from Cython:

- ``scipy.optimize.bisect``
- ``scipy.optimize.ridder``
- ``scipy.optimize.brenth``
- ``scipy.optimize.brentq``

The Cython Optimize Zeros API could be used to "vectorize" a root finder given
a single callback function and an array of input conditions using a loop in C
or with Cython ``prange``. Import the module into Cython as follows::

    from scipy.optimize cimport cython_optimize


Callback Signature
------------------
There is currently only one callback signature that can be used with the root
finders accessible from Cython.

    double (*callback_type)(double, void*)


``cython_optimize.zeros_struct``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a user defined ``struct``
with any extra parameters as the 2nd argument. ::

Import the module into Cython as follows::

    from scipy.optimize.cython_optimize cimport zeros_struct


Examples
--------
Usage of ``scipy.optimize.cython_optimize`` requires the use of Cython to write
callbacks that are compiled into C.

These are the basic steps to use ``scipy.optimize.cython_optimize``:

1. Create a Cython ``.pyx`` file
2. Import the ``zeros_struct`` module
3. Write the callback function, and call the selected root-finder passing the
   callback, any extra arguments, and the other solver parameters ::

       import math
       from scipy.optimize.cython_optimize cimport zeros_struct

       ARGS = {'C0': 1.0, 'C1': 0.7}
       XTOL, RTOL, MITR = 1e-3, 1e-3, 10

       ctypedef struct test_params:
           double C0
           double C1


       cdef double f(double x, void *args):
           cdef test_params *myargs = <test_params *> args
           return myargs.C0 - math.exp(-(x - myargs.C1))


       cdef double brentq_wrapper_example(dict args, double xa, double xb,
                                          double xtol, double rtol, int mitr):
           # Cython automatically casts dictionary to struct
           cdef test_params myargs = args  
           return zeros_struct.brentq(
               f, xa, xb, <test_params *> &myargs, xtol, rtol, mitr, NULL)


       def brentq_example(args=ARGS, xa=0.5, xb=1.0, xtol=XTOL, rtol=RTOL, mitr=MITR):
           return brentq_wrapper_example(args, xa, xb, xtol, rtol, mitr)


       x = brentq_example()
       # 0.6999942848231314

4. If you want to call your function from Python, create a wrapper
5. Create a Cython ``.pxd`` file if you need to export any Cython functions


Full Output
-----------
The  functions in ``scipy.optimize.cython_optimize`` can also copy the full
output from the solver to a C ``struct`` that is passed as its last argument.
If you don't want the full output just pass ``NULL``. The full output
``struct`` should contain the following:

- ``funcalls``: number of function calls
- ``iterations``: number of iterations
- ``error_num``: error number

The error number returned is -1 for a sign error, and -2 for a convergence
error. All other values mean the solver converged. Note that the full output
``struct`` must be cast to ``scipy_zeros_parameters``. ::

    from scipy.optimize.cython_optimize cimport zeros_struct
    from scipy.optimize.cython_optimize.tests.examples.zeros_struct_examples \
        cimport f_solarcell

    ARGS = {'dark_current': 1e-09, 'series_resistance': 0.004,
            'shunt_resistance': 10.0, 'thermal_voltage': 0.27456}
    XTOL, RTOL, MITR = 0.001, 0.001, 10

    ctypedef struct test_params:
        double voltage
        double light_current
        double dark_current
        double series_resistance
        double shunt_resistance
        double thermal_voltage

    ctypedef struct scipy_brent_full_output:
        int funcalls
        int iterations
        int error_num
        double root

    # cython brentq solver with full output
    cdef scipy_brent_full_output solarcell_brent_full_output(
            dict args, double xa, double xb, double xtol, double rtol,
            int mitr):
        cdef test_params myargs = args
        cdef scipy_brent_full_output full_output
        full_output.root = zeros_struct.brentq(
            f_solarcell, xa, xb, <test_params *> &myargs, xtol, rtol, mitr,
            <zeros_struct.scipy_zeros_parameters *> &full_output)
        return full_output


    def test_brent_full_output(args=ARGS, xa=0.0, xb=6.0, xtol=XTOL, rtol=RTOL,
                               mitr=MITR):
        '''test brent with full output'''
        return solarcell_brent_full_output(args, xa, xb, xtol, rtol, mitr)

    result = test_brent_full_output()
    # {'funcalls': 4,
    #  'iterations': 3,
    #  'error_num': 281790720,
    #  'root': 5.255231961257658}
"""
