"""
Cython Optimize Zeros API
=========================
This package contains a Cython API for the root finding scalar functions in
:py:mod:`scipy.optimize` to obtain efficient pure C implementations. Faster
execution is not guaranteed. The relative performance of ``optimize`` and
``cython_optimize`` depends on the application and ranges from slower to 10X
faster. The Cython Optimize Zeros API could be used for native C looping or
with Cython ``prange`` for efficient execution of large arrays of inputs with
the same callbacks.

Callback Signatures
-------------------
The Cython Optimize Zeros API provides three callback signagures for the root
finders, two which are safe without the global interpreter lock (GIL) and one
that uses Python tuples.

- :py:mod:`scipy.optimize.cython_optimize.zeros_struct`::

        ``double (*callback_type)(double, void*)``

  This callback takes a double with the scalar independent variable as the 1st
  argument and a user defined ``struct`` with any extra parameters as the 2nd.

- :py:mod:`scipy.optimize.cython_optimize.zeros_array`::

        ``double (*callback_type_array)(int, double*)``

  This callback takes an integer with the number of extra parameters as the 1st
  argument and an array of doubles with any extra parameters as the 2nd. Even
  if the integer is unused in your callback, it must still be in the signature,
  because an internal wrapper will use it to prepend a double with the scalar
  independent variable to the array. In your callback, the independent variable
  must be the first element in the array, followed by the extra parameters.

- :py:mod:`scipy.optimize.cython_optimize.zeros_tuple`::

        ``double (*callback_type_tuple)(double, tuple)``

  This callback takes a double with the scalar independent variable as the 1st
  argument and a Python tuple with any extra parameters as the 2nd. Therefore
  this signature is not safe to use with ``nogil``. Also the maximum number of
  extra parameters is hard-coded as ``MAXARGS = 10``.

Available Functions
-------------------
There are seven available functions which are all implemented in pure C, and
exposed in modules matching the expected callback signature.

* ``scipy.optimize.cython_optimize.zeros_struct`` exposes root finders expecting
  extra arguments in a native C ``struct``. See examples for usage.

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.newton`::

            double newton(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.secant`::

            double secant(callback_type func, double x0, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.halley`::

            double halley(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, callback_type fprime2, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.bisect`::

            double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.ridder`::

            double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.brenth`::

            double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_struct.brentq`::

            double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

* ``scipy.optimize.cython_optimize.zeros_array`` exposes root finders expecting
  extra arguments in a native C array of doubles. See examples for usage.

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.newton`::

            double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.secant`::

            double secant(callback_type_array func, double x0, int n, double* args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.halley`::

            double halley(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.bisect`::

            double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.ridder`::

            double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.brenth`::

            double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_array.brentq`::

            double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

* ``scipy.optimize.cython_optimize.zeros_tuple`` exposes root finders expecting extra
  arguments in a Python tuple.

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.newton`::

            double newton(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.secant`::

            double secant(callback_type_tuple func, double x0, tuple args, double tol, int maxiter, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.halley`::

            double halley(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, callback_type_tuple fprime2, scipy_newton_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.bisect`::

            double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.ridder`::

            double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.brenth`::

            double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

    - :py:func:`scipy.optimize.cython_optimize.zeros_tuple.brentq`::

            double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

Examples
--------
Usage of the Cython Optimize Zeros API requires the use of Cython to write
callbacks that are compiled into native C. See the ``cython_optimize`` Jupyter
notebook in the ``cython_optimize/Examples`` folder for basic usage examples.
There are also several realistic examples that solve the single-diode model of
a solar cell in the ``cython_optimize/Examples`` folder. These examples are
also used to test the different combinations of callbacks and root-finders.

These are the basic steps to use ``cython_optimize``:

1. Select a callback signature, *eg* ``scipy.optimize.cython_optimize.zeros_struct``
2. Select the root finder, *eg* :py:func:`~scipy.optimize.cython_optimize.zeros_struct.newton`
3. Create a Cython ``.pyx`` file that imports the zeros module with the
   selected callback signature, create the callback(s) using the selected
   signature, and call the selected root-finder passing the callback(s) and any
   extra arguments, and other parmetres paramters
4. If you want to call your function from Python, create a wrapper
5. Create a Cython ``.pxd`` file if you need to export any Cython functions
"""

from __future__ import division, print_function, absolute_import

__all__ = ['zeros_tuple', 'zeros_struct', 'zeros_array']

from . import zeros_tuple, zeros_struct, zeros_array

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
