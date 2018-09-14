"""
Cython Optimize Zeros API
=========================
This package contains a Cython API for some of the root finding scalar
functions in :py:mod:`scipy.optimize` to obtain efficient pure C
implementations, although faster execution is not guaranteed. For example the
Cython Optimize Zeros API could be used for native C looping or with Cython
``prange`` for efficient execution of large arrays of inputs with the same
callbacks. The Cython Optimize Zeros API can be imported into Cython code using
the following::

    from scipy.optimize cimport cython_optimize

Callback Signatures
-------------------
The Cython Optimize Zeros API provides three callback signagures for the root
finders, two which are safe without the global interpreter lock (GIL) and one
that uses Python tuples. The functions are grouped by callback signatures into
Cython modules that can be imported using ``cimport``.

.. py:module:: scipy.optimize.cython_optimize.zeros_struct

``cython_optimize.zeros_struct``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a user defined ``struct``
with any extra parameters as the 2nd. ::

    double (*callback_type)(double, void*)

This module can be imported using the following::

    from scipy.optimize.cython_optimize cimport zeros_struct

.. py:module:: scipy.optimize.cython_optimize.zeros_array

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

    double (*callback_type_array)(int, double*)

This module can be imported using the following::

    from scipy.optimize.cython_optimize cimport zeros_array

.. py:module:: scipy.optimize.cython_optimize.zeros_tuple

``cython_optimize.zeros_tuple``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The root-finders in this module use a callback that takes a double with the
scalar independent variable as the 1st argument and a Python tuple with any
extra parameters as the 2nd. Therefore this signature is *not* safe to use with
``nogil``. ::

    double (*callback_type_tuple)(double, tuple)

This module can be imported using the following::

    from scipy.optimize.cython_optimize cimport zeros_tuple

Available Functions
-------------------
There are four available functions which are all implemented in pure C and
exposed in modules matching the expected callback signature. These root-finding
functions correspond to the functions available in ``scipy.optimize.zeros``.

- :py:func:`~scipy.optimize.bisect`
- :py:func:`~scipy.optimize.ridder`
- :py:func:`~scipy.optimize.brenth`
- :py:func:`~scipy.optimize.brentq`

:py:mod:`~scipy.optimize.cython_optimize.zeros_struct`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module exposes root finders expecting extra arguments in a native C
``struct``. See examples for usage.

- ``cython_optimize.zeros_struct.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.bisect(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f``, that
     takes a double and a native C ``struct`` ::

         double bisect(callback_type f, double xa, double xb, void* args,
                       double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_struct.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.ridder(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar to :py:func:`~scipy.optimize.ridder` with a callback, ``f``, that
     takes a double and a native C ``struct`` ::

         double ridder(callback_type f, double xa, double xb, void* args,
                       double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_struct.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.brenth(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar to :py:func:`~scipy.optimize.brenth` with a callback, ``f``, that
     takes a double and a native C ``struct`` ::

         double brenth(callback_type f, double xa, double xb, void* args,
                       double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_struct.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_struct.brentq(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similar to :py:func:`~scipy.optimize.brentq` with a callback, ``f``, that
     takes a double and a native C ``struct`` ::

         double brentq(callback_type f, double xa, double xb, void* args,
                       double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param void* args: a ``struct`` of extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

:py:mod:`~scipy.optimize.cython_optimize.zeros_array`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exposes root finders expecting extra arguments in a native C array of doubles.
See examples for usage.

- ``cython_optimize.zeros_array.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.bisect(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f``, that
     takes an integer and an array of doubles ::

        double bisect(callback_type_array f, double xa, double xb, int n, double* args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_array.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.ridder(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.ridder` with a callback, ``f``, that
     takes an integer and an array of doubles ::

        double ridder(callback_type_array f, double xa, double xb, int n, double* args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_array.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.brenth(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brenth` with a callback, ``f``, that
     takes an integer and an array of doubles ::

        double brenth(callback_type_array f, double xa, double xb, int n, double* args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_array.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_array.brentq(f, xa, xb, n, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brentq` with a callback, ``f``, that
     takes an integer and an array of doubles ::

        double brentq(callback_type_array f, double xa, double xb, int n, double* args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_array f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param int n: the number of extra parameters
     :param double* args: an array of extra parameters, max size of 10
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

:py:mod:`~scipy.optimize.cython_optimize.zeros_tuple`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exposes root finders expecting extra arguments in a Python tuple.

- ``cython_optimize.zeros_tuple.bisect``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.bisect(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.bisect` with a callback, ``f``, that
     takes a double and a Python tuple ::

        double bisect(callback_type_tuple f, double xa, double xb, tuple args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_tuple.ridder``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.ridder(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.ridder` with a callback, ``f``, that
     takes a double and a Python tuple ::

        double ridder(callback_type_tuple f, double xa, double xb, tuple args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_tuple.brenth``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.brenth(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brenth` with a callback, ``f``, that
     takes a double and a Python tuple ::

        double brenth(callback_type_tuple f, double xa, double xb, tuple args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

- ``cython_optimize.zeros_tuple.brentq``

  .. py:function:: scipy.optimize.cython_optimize.zeros_tuple.brentq(f, xa, xb, args, xtol, rtol, iter, full_output)

     Similiar to :py:func:`~scipy.optimize.brentq` with a callback, ``f``, that
     takes a double and a Python tuple ::

        double brentq(callback_type_tuple f, double xa, double xb, tuple args,
                      double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

     :param callback_type_tuple f: the callback function
     :param double xa: lower bound
     :param double xb: upper bound
     :param tuple args: extra parameters
     :param double xtol: absolute tolerance
     :param double rtol: relative tolerance
     :param int iter: maximum number of iterations
     :param scipy_zeros_parameters* full_output: optional output ``struct``, use ``NULL`` if not desired
     :return: double

Examples
--------
Usage of :py:mod:`~scipy.optimize.cython_optimize` requires the use of Cython
to write callbacks that are compiled into native C. There are several usage examples that solve the single-diode model of
a solar cell in the ``cython_optimize/Examples`` folder. These examples are
also used to test the different combinations of callbacks and root-finders.

These are the basic steps to use :py:mod:`~scipy.optimize.cython_optimize`:

1. Select a callback signature, for example:
   :py:mod:`~scipy.optimize.cython_optimize.zeros_struct`
2. Select the root finder, for example:
   :py:func:`~scipy.optimize.cython_optimize.zeros_struct.brentq`
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
The  functions in :py:mod:`~scipy.optimize.cython_optimize` can also copy the
full output from the solver to a C ``struct`` that is passed as its last
argument. If you don't want the full output just pass ``NULL``. The full output
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
