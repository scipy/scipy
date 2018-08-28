"""
Cython Optimize Zeros API
=========================
This package contains a Cython API for the root finding scalar functions in
:py:mod:`scipy.optimize` to obtain efficient pure C implementations, although
faster execution is not guaranteed. For example the Cython Optimize Zeros API
could be used for native C looping or with Cython ``prange`` for efficient
execution of large arrays of inputs with the same callbacks. The Cython
Optimize Zeros API can be imported into Cython code using the following::

    cimport scipy.optimize.cython_optimize as coz

Callback Signatures
-------------------
The Cython Optimize Zeros API provides three callback signagures for the root
finders, two which are safe without the global interpreter lock (GIL) and one
that uses Python tuples. The functions are grouped by callback signatures into
Cython modules that can be imported using ``cimport``.

.. py:module:: scipy.optimize.cython_optimize.zeros_struct

``scipy.optimize.cython_optimize.zeros_struct``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module contains a callback takes a double with the scalar independent
variable as the 1st argument and a user defined ``struct`` with any extra
parameters as the 2nd. ::

    double (*callback_type)(double, void*)

.. py:module:: scipy.optimize.cython_optimize.zeros_array

``scipy.optimize.cython_optimize.zeros_array``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module contains a callback takes an integer with the number of extra
parameters as the 1st argument and an array of doubles with any extra
parameters as the 2nd. Even if the integer is unused in your callback, it must
still be in the signature, because an internal wrapper will use it to prepend a
double with the scalar independent variable to the array. In your callback, the
independent variable must be the first element in the array, followed by the
extra parameters.  Also the maximum number of extra parameters is hard-coded as
``MAXARGS = 10``. ::

    double (*callback_type_array)(int, double*)

.. py:module:: scipy.optimize.cython_optimize.zeros_tuple

``scipy.optimize.cython_optimize.zeros_tuple``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module contains a callback takes a double with the scalar independent
variable as the 1st argument and a Python tuple with any extra parameters as
the 2nd. Therefore this signature is not safe to use with ``nogil``. ::

    double (*callback_type_tuple)(double, tuple)

Available Functions
-------------------
There are seven available functions which are all implemented in pure C and
exposed in modules matching the expected callback signature. These root-finding
functions correspond to the functions available in ``scipy.optimize.zeros``.

:py:mod:`~scipy.optimize.cython_optimize.zeros_struct`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module exposes root finders expecting extra arguments in a native C
``struct``. See examples for usage.

- ``scipy.optimize.cython_optimize.zeros_struct.newton``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.newton(func, x0, fprime, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with two callbacks, ``func``
     and ``fprime``, that take a double and a native C ``struct`` ::

         double newton(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type func: the callback function
     :param double x0: the initial guess
     :param callback_type fprime: the first derivative of the callback
     :param void* args: a ``struct`` of extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_struct.secant``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.secant(func, x0, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with one callback, ``func``
     that takes a double and a native C ``struct`` ::

         double secant(callback_type func, double x0, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type func: the callback function
     :param double x0: the initial guess
     :param void* args: a ``struct`` of extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_struct.halley``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.halley(func, x0, fprime, args, tol, maxiter, fprime2, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with three callbacks,
     ``func``, ``fprime``, and ``fprime2`` that each takes double and a native
     C ``struct`` ::

         double halley(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, callback_type fprime2, scipy_newton_parameters *full_output)

     :param callback_type func: the callback function
     :param double x0: the initial guess
     :param callback_type fprime: the first derivative of the callback
     :param void* args: a ``struct`` of extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param callback_type fprime2: the second derivative of the callback
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_struct.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.bisect(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f`` that
     takes a double and a native C ``struct`` ::

         double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- `scipy.optimize.cython_optimize.zeros_struct.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.ridder(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar :py:func:`~scipy.optimize.ridder` with a callback, ``f`` that
     takes a double and a native C ``struct`` ::

         double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_struct.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.brenth(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar to :py:func:`~scipy.optimize.brenth` with a callback, ``f`` that
     takes a double and a native C ``struct`` ::

         double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_struct.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.brentq(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar to :py:func:`~scipy.optimize.brentq` with a callback, ``f`` that
     takes a double and a native C ``struct`` ::

         double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

:py:mod:`~scipy.optimize.cython_optimize.zeros_array`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exposes root finders expecting extra arguments in a native C array of doubles.
See examples for usage.

- ``scipy.optimize.cython_optimize.zeros_array.newton``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.newton(func, x0, fprime, n, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with two callbacks, ``func``
     and ``fprime``, that take an integer and an array of doubles ::

         double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type_array func: the callback function
     :param double x0: the initial guess
     :param callback_type_array fprime: the first derivative of the callback
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_array.secant``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.secant(func, x0, n, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with one callback, ``func``
     that takes an integer and an array of doubles ::

        double secant(callback_type_array func, double x0, int n, double* args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type_array func: the callback function
     :param double x0: the initial guess
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_array.halley``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.halley(func, x0, fprime, n, args, tol, maxiter, fprime2, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with three callbacks,
     ``func``, ``fprime``, and ``fprime2`` that each takes an integer and an
     array of doubles ::

        double halley(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2, scipy_newton_parameters *full_output)

     :param callback_type_array func: the callback function
     :param double x0: the initial guess
     :param callback_type_array fprime: the first derivative of the callback
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param callback_type_array fprime2: the second derivative of the callback
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_array.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.bisect(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f`` that
     takes an integer and an array of doubles ::

        double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_array.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.ridder(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.ridder` with a callback, ``f`` that
     takes an integer and an array of doubles ::

        double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_array.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.brenth(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brenth` with a callback, ``f`` that
     takes an integer and an array of doubles ::

        double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_array.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.brentq(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brentq` with a callback, ``f`` that
     takes an integer and an array of doubles ::

        double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

:py:mod:`~scipy.optimize.cython_optimize.zeros_tuple`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exposes root finders expecting extra arguments in a Python tuple.

- ``scipy.optimize.cython_optimize.zeros_tuple.newton``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.newton(func, x0, fprime, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with two callbacks, ``func``
     and ``fprime``, that take a double and a Python tuple ::

        double newton(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type_tuple func: the callback function
     :param double x0: the initial guess
     :param callback_type_tuple fprime: the first derivative of the callback
     :param tuple args: extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_tuple.secant``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.secant(func, x0, args, tol, maxiter, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with one callback, ``func``
     that takes a double and a Python tuple ::

        double secant(callback_type_tuple func, double x0, tuple args, double tol, int maxiter, scipy_newton_parameters *full_output)

     :param callback_type_tuple func: the callback function
     :param double x0: the initial guess
     :param tuple args: extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_tuple.halley``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.halley(func, x0, fprime, args, tol, maxiter, fprine2, full_ouptput)

     Similar to :py:func:`~scipy.optimize.newton` with three callbacks,
     ``func``, ``fprime``, and ``fprime2`` that each takes a double and a
     Python tuple ::

        double halley(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, callback_type_tuple fprime2, scipy_newton_parameters *full_output)

     :param callback_type_tuple func: the callback function
     :param double x0: the initial guess
     :param callback_type_tuple fprime: the first derivative of the callback
     :param tuple args: extra parameters
     :param double tol: exit tolerance
     :param int maxiter: maximum number of iterations
     :param callback_type_tuple fprime2: the second derivative of the callback
     :param scipy_newton_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :returns: double

- ``scipy.optimize.cython_optimize.zeros_tuple.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.bisect(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f`` that
     takes a double and a Python tuple ::

        double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_tuple.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.ridder(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.ridder` with a callback, ``f`` that
     takes a double and a Python tuple ::

        double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_tuple.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.brenth(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brenth` with a callback, ``f`` that
     takes a double and a Python tuple ::

        double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

- ``scipy.optimize.cython_optimize.zeros_tuple.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.brentq(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brentq` with a callback, ``f`` that
     takes a double and a Python tuple ::

        double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, set to ``NULL`` if full output is not desired
     :return: double

Examples
--------
Usage of the Cython Optimize Zeros API requires the use of Cython to write
callbacks that are compiled into native C. See the ``cython_optimize`` Jupyter
notebook in the ``cython_optimize/Examples`` folder for basic usage examples.
There are also several realistic examples that solve the single-diode model of
a solar cell in the ``cython_optimize/Examples`` folder. These examples are
also used to test the different combinations of callbacks and root-finders.

These are the basic steps to use ``cython_optimize``:

1. Select a callback signature, *eg* :py:mod:`~scipy.optimize.cython_optimize.zeros_struct`
2. Select the root finder, *eg* :py:func:`~scipy.optimize.cython_optimize.zeros_struct.newton`
3. Create a Cython ``.pyx`` file that imports the zeros module with the
   selected callback signature, create the callback(s) using the selected
   signature, and call the selected root-finder passing the callback(s) and any
   extra arguments, and other parmetres paramters ::

       import math
       cimport scipy.optimize.cython_optimize as coz

       XTOL, RTOL, MITR = 1e-3, 1e-3, 10

       ctypedef struct test_params:
           double C0
           double C1


       cdef double f(double x, void *args):
           cdef test_params *myargs = <test_params *> args
           return myargs.C0 - math.exp(-(x - myargs.C1))


       cdef test_params myargs
       myargs.C0 = 1.0
       myargs.C1 = 0.7

       x = coz.zeros_struct.brentq(f, 0.5, 1, <test_params *> &myargs, XTOL, RTOL, MITR, NULL))
       # 0.6999942848231314

4. If you want to call your function from Python, create a wrapper
5. Create a Cython ``.pxd`` file if you need to export any Cython functions

Full Output
-----------
The  Cython Optimize Zeros API can also copy the full output from the solver.
Create a  C ``struct`` and pass it as the last argument. If you don't want the
full output just pass ``NULL`` as the last argument.

- number of function calls: ``funcalls``
- number of iterations: ``iterations``
- error number: ``error_num``

The ``struct`` must be the cast to ``scipy_newton_parameters`` for the three
newton functions and ``scipy_zeros_parameters`` for the four zeros functions. ::

    from scipy.optimize.cython_optimize cimport zeros_tuple
    from scipy.optimize.cython_optimize.examples.zeros_tuple_examples cimport f_solarcell, fprime, fprime2

    ARGS = (5.25, 6.0, 1e-09, 0.004, 10.0, 0.27456)
    XTOL, RTOL, MITR = 0.001, 0.001, 10

    DEF SIGNERR = -1
    DEF CONVERR = -2

    ctypedef struct scipy_brent_full_output:
        int funcalls
        int iterations
        int error_num
        double root

    # cython brentq solver with full output
    cdef scipy_brent_full_output solarcell_brent_full_output(tuple args, double xa, double xb, double xtol, double rtol, int mitr):
        cdef scipy_brent_full_output full_output
        full_output.root = zeros_tuple.brentq(f_solarcell, xa, xb, args, xtol, rtol, mitr, <zeros_tuple.scipy_zeros_parameters *> &full_output)
        return full_output


    def test_brent_full_output(args=ARGS, xa=0.0, xb=6.0, xtol=XTOL, rtol=RTOL, mitr=MITR):
        '''test newton with full output'''
        return solarcell_brent_full_output(args, xa, xb, xtol, rtol, mitr)

    result = test_brent_full_output()
    # {'funcalls': 4,
    #  'iterations': 3,
    #  'error_num': 281790720,
    #  'root': 5.255231961257658}
"""

from __future__ import division, print_function, absolute_import

__all__ = ['zeros_tuple', 'zeros_struct', 'zeros_array']

from . import zeros_tuple, zeros_struct, zeros_array

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
