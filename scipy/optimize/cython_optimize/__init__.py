"""
Cython Optimize Zeros API
=========================
The following four zeros functions can be accessed directly from Cython:

- ``bisect``
- ``ridder``
- ``brenth``
- ``brentq``

Import the module into Cython as follows::

    from scipy.optimize cimport cython_optimize

The Cython Optimize Zeros API could be used to "vectorize" a root finder given
a single callback function and an array of input conditions using a for-loop or
``prange`` in Cython.


Callback Signature
------------------
The zeros functions in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a user defined ``struct``
with any extra parameters as the 2nd argument. ::

    double (*callback_type)(double, void*)

Import the module into Cython as follows::

    from scipy.optimize.cython_optimize cimport zeros


Examples
--------
Usage of ``scipy.optimize.cython_optimize`` requires Cython to write callbacks
that are compiled into C. For more information on compiling Cython see the
`Cython Documentation <http://docs.cython.org/en/latest/index.html>`_.

These are the basic steps:

1. Create a Cython ``.pyx`` file
2. Import the ``cython_optimize.zeros`` module
3. Write the callback function, and call the selected zeros function passing
   the callback, any extra arguments, and the other solver parameters ::

       import math
       from scipy.optimize.cython_optimize cimport zeros

       ARGS = {'C0': 1.0, 'C1': 0.7}  # a dictionary of extra arguments
       XLO, XHI = 0.5, 1.0  # lower and upper search boundaries
       XTOL, RTOL, MITR = 1e-3, 1e-3, 10  # other solver parameters

       # user defined struct for extra parameters
       ctypedef struct test_params:
           double C0
           double C1


       # user defined callback
       cdef double f(double x, void *args):
           cdef test_params *myargs = <test_params *> args
           return myargs.C0 - math.exp(-(x - myargs.C1))


       # Cython wrapper function
       cdef double brentq_wrapper_example(dict args, double xa, double xb,
                                          double xtol, double rtol, int mitr):
           # Cython automatically casts dictionary to struct
           cdef test_params myargs = args  
           return zeros.brentq(
               f, xa, xb, <test_params *> &myargs, xtol, rtol, mitr, NULL)


       # Python function
       def brentq_example(args=ARGS, xa=XLO, xb=XHI, xtol=XTOL, rtol=RTOL,
                          mitr=MITR):
           '''Calls Cython wrapper from Python.'''
           return brentq_wrapper_example(args, xa, xb, xtol, rtol, mitr)


       x = brentq_example()
       # 0.6999942848231314

4. If you want to call your function from Python, create a Cython wrapper, and
   a Python function that calls the wrapper.
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

An error number of -1 means a sign error, -2 means a convergence error, and 0
means the solver converged. Note that the full output ``struct`` must be cast
to ``scipy_zeros_parameters``. Continuing from the previous example::


    # user defined full output structure with required fields
    ctypedef struct scipy_brent_full_output:
        int funcalls
        int iterations
        int error_num
        double root


    # cython brentq solver with full output
    cdef scipy_brent_full_output brentq_full_output_wrapper_example(
            dict args, double xa, double xb, double xtol, double rtol,
            int mitr):
        cdef test_params myargs = args
        cdef scipy_brent_full_output full_output  # used instead of NULL
        # put result into full_output
        full_output.root = zeros.brentq(
            f, xa, xb, <test_params *> &myargs, xtol, rtol, mitr,
            <zeros.scipy_zeros_parameters *> &full_output)
        return full_output


    # Python function
    def brent_full_output_example(args=ARGS, xa=XLO, xb=XHI, xtol=XTOL,
                                  rtol=RTOL, mitr=MITR):
        '''Returns full output'''
        return brentq_full_output_wrapper_example(args, xa, xb, xtol, rtol,
                                                  mitr)

    result = brent_full_output_example()
    # {'error_num': 0,
    #  'funcalls': 6,
    #  'iterations': 5,
    #  'root': 0.6999942848231314}
"""
