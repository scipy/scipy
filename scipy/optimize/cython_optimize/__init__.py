"""
Cython Optimize Zeros API
=========================
This package contains a Cython API for the root finding :ref:`scalar-functions`
in :py:mod:`scipy.optimize` to obtain performant and efficient pure C
implementations. Faster execution is not guaranteed, but you might get 5-10X
improvement versus pure Python. The Cython Optimize Zeros API  could be used
for native C looping or with Cython `prange` for efficient execution of large
arrays of inputs with the same callbacks.

Callback Signatures
-------------------
The Cython Optimize Zeros API provides three callback signagures for the root
finders, two which are safe without the global interpreter lock (GIL) and one
that uses Python tuples.

- ``double (*callback_type)(double, void*)``
- ``double (*callback_type_array)(int, double*)``
- ``double (*callback_type_tuple)(double, tuple)``

Available Functions
-------------------
There are seven available functions which are all implemented in pure C, and
exposed in modules matching the expected callback signature.

* ``scipy.optimize.cython_optimize.zeros_struct`` exposes root finders expecting
  extra arguments in a native C ``struct``. See examples for usage.

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.newton`::

            cdef double newton(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.secant`::

            double secant(callback_type func, double x0, void* args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.halley`::

            double halley(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, callback_type fprime2)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.bisect`::

            double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.ridder`::

            double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.brenth`::

            double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_struct.brentq`::

            double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

* ``scipy.optimize.cython_optimize.zeros_array`` exposes root finders expecting
  extra arguments in a native C array of doubles. See examples for usage.

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.newton`::

            double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.secant`::

            double secant(callback_type_array func, double x0, int n, double* args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.halley`::

            double halley(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.bisect`::

            double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.ridder`::

            double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.brenth`::

            double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros_array.brentq`::

            double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

* ``scipy.optimize.cython_optimize.zeros`` exposes root finders expecting extra
  arguments in a Python tuple.

    - :py:func:`~scipy.optimize.cython_optimize.zeros.newton`::

            double newton(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.secant`::

            double secant(callback_type_tuple func, double x0, tuple args, double tol, int maxiter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.halley`::

            double halley(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, callback_type_tuple fprime2)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.bisect`::

            double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.ridder`::

            double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.brenth`::

            double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

    - :py:func:`~scipy.optimize.cython_optimize.zeros.brentq`::

            double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

Examples
--------
Usage of the Cython Optimize Zeros API requires the use of Cython to write
callbacks that are compiled into native C. There are several usage examples
with different combinations of callbacks and root-finders in the
``cython_optimize/Examples`` folder. There are five steps:

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

__all__ = ['zeros', 'zeros_struct', 'zeros_array']

from . import zeros, zeros_struct, zeros_array

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
