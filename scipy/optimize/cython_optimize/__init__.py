"""
Cython Optimize
===============
This package contains an experimental Cython API for ``scipy.optimize``.
This is a WIP. Only ``zeros.newton`` and ``zeros.bisect`` are currently
supported.

Proposal
--------
The Cython API will not use any Python objects so that it can be used GIL free
if desired. The callback will use the basic
``double (*callback_type)(double,void*)`` type definition.

Background
----------
Currently there are 5 methods in ``zeros``:

    * `newton`
    * `bisect`
    * `ridder`
    * `brentq`
    * `brenth`

Only ``newton`` is pure Python, and is contained in the ``zeros.py`` module
itself. The other four are written in native C with no Python bindings. This is
advantageous because they can be used GIL free if the Cython API is also Python
free.

``zeros.c``
~~~~~~~~~~~
The four native C functions are currently exposed using the Python C-API in
``zeros.c``, which is exported as the private Python module ``_zeros`` and
imported into ``zeros.py``. The Python methods in ``zeros.py`` are mapped to
private methods in ``zeros.c`` called ``_bisect``, ``_ridder``, ``_brentq``, and
``_brenth``. These methods all execute the function ``call_solver``, which is
_not_ exported. As the first argument they all pass a ``solver_type``, which is
a reference to the native C function. The other argument is a tuple of all of
the arguments from the original caller in ``zeros.py``.

``call_solver``
~~~~~~~~~~~~~~~
This hidden method performs the redundant steps of parsing the arguments tuple,
error checking, allocating a new tuple to include the independent variable, and
incrementing the reference count of each stolen item from the original tuple of
additional callback function parameters.

The parsed argument tuple is the same for all four methods:

* ``f`` - Python callback function that takes a tuple
* ``a`` - lower limit of the independent variable
* ``b`` - upper limit of the independent variable
* ``xtol`` - absolute tolerance
* ``rtol`` - relative tolerance
* ``iter`` - maximum iterations
* ``xargs`` - Python tuple of extra callback arguments
* ``fulloutput`` - flag to use ``RootResults`` object which returns iterations,
function calls, and converged flag
* ``disp`` - raise ``RuntimeError`` if didn't converge

``scipy_zeros_parameters``
~~~~~~~~~~~~~~~~~~~~~~~~~~
The Python module ``zeros.c`` mimics the ``default_parameters`` structure found
in the ``zeros.h`` header, but adds on 3 python objects, and calls the new
structure, ``scipy_zeros_parameters``:

* Python callback function passed to the ``zeros.py`` method
* Python tuple containing the extra callback arguments.
* a native C jump buffer called ``env`` used to quit execution and return
``NULL`` if the Python callback returns ``NULL``.

The arguments are copied to a new temporary Python tuple one item larger for the
independent variable, and the copied references are incremented, so that when
the temporary tuple is decremeneted, the arguments are still referenced. Then
the ``scipy_zeros_parameters`` structure fields are set to the Python callback
``f`` and the temporary ``fargs`` tuple. Then if ``setjmp`` returns 0, _ie_: it
was called directly and not the result of ``longjmp``, the native C solver
function is called, with a native C callback function called
``scipy_zeros_functions_func`` as the first argument.

``scipy_zeros_functions_func``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a function wrapper. It mimics ``callback_type`` from the ``zeros.h``
header, which takes a double and a pointer to a ``default_parameters``
structure. When the wrapper executes, it gets the Python callback function and
the the arguments tuple from the structure, sets the first item in the tuple to
the double, passes it the callback, and executes it. If the callback fails, then
it executes a long jump to the set jump point, which terminates ``call_solver``,
decrements the ``fargs`` temporary tuple, and returns ``NULL``. Otherwise, the
wrapper function copies the return value to a native double, decrements the
Python return value object, and returns the native double value to the solver.

The solver continues iterating and ultimately returns the zero to the Python
module, which either copies it to the ``RootResults`` object, or copies it to a Python
float, and returns.
"""

from __future__ import division, print_function, absolute_import

__all__ = ['zeros', 'zeros_struct']

from . import zeros, zeros_struct

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
