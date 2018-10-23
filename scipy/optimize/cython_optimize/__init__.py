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

Callback Signatures
-------------------
Three different callback signatures can be used with the root finders. The functions are grouped by callback signatures into
Cython modules that can be imported using ``cimport``.

    double (*callback_type)(double, void*)
    double (*callback_type_array)(int, double*)
    double (*callback_type_tuple)(double, tuple)

Import the module with the desired callback type.

``cython_optimize.zeros_struct``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a user defined ``struct``
with any extra parameters as the 2nd. ::

Import the module into Cython as follows::

    from scipy.optimize.cython_optimize cimport zeros_struct

``cython_optimize.zeros_array``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes an integer with the
number of extra parameters as the 1st argument and an array of doubles with any
extra parameters as the 2nd. Even if the integer is unused in your callback, it
must still be in the signature, because an internal wrapper will use it to
prepend the scalar independent variable to the array, which is then passed to
your callback. Therefore, in your callback, the independent variable must be
the first element in the array, followed by the extra parameters.  Also the
maximum number of extra parameters is hard-coded as ``MAXARGS = 10``. ::

Import the module into Cython as follows::

    from scipy.optimize.cython_optimize cimport zeros_array

``cython_optimize.zeros_tuple``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a Python tuple with any
extra parameters as the 2nd. Therefore this signature is *not* safe to use with
``nogil``. ::

Import the module into Cython as follows::

    from scipy.optimize.cython_optimize cimport zeros_tuple

Examples
--------
Usage of ``scipy.optimize.cython_optimize`` requires the use of Cython to write
callbacks that are compiled into C.

These are the basic steps to use ``scipy.optimize.cython_optimize``:

1. Select a callback signature, for example:
   ``scipy.optimize.cython_optimize.zeros_struct``
2. Select the root finder, for example:
   ``scipy.optimize.cython_optimize.zeros_struct.brentq``
3. Create a Cython ``.pyx`` file that imports the zeros module with the
   selected callback signature, create the callback(s) using the selected
   signature, and call the selected root-finder passing the callback(s), any
   extra arguments, and the other solver parameters ::

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
           cdef test_params myargs = args  # Cython automatically casts dictionary to struct
           return zeros_struct.brentq(f, xa, xb, <test_params *> &myargs, xtol, rtol, mitr, NULL)


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

    from scipy.optimize.cython_optimize cimport zeros_tuple
    from scipy.optimize.cython_optimize.examples.zeros_tuple_examples cimport f_solarcell

    ARGS = (5.25, 6.0, 1e-09, 0.004, 10.0, 0.27456)
    XTOL, RTOL, MITR = 0.001, 0.001, 10

    ctypedef struct scipy_brent_full_output:
        int funcalls
        int iterations
        int error_num
        double root

    # cython brentq solver with full output
    cdef scipy_brent_full_output solarcell_brent_full_output(
            tuple args, double xa, double xb, double xtol, double rtol,
            int mitr):
        cdef scipy_brent_full_output full_output
        full_output.root = zeros_tuple.brentq(
            f_solarcell, xa, xb, args, xtol, rtol, mitr,
            <zeros_tuple.scipy_zeros_parameters *> &full_output)
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
