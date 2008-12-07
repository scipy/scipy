General information
===================

Examples in this tutorial
-------------------------

Throughout this tutorial it is assumed that the user
has imported all of the names defined in the SciPy top-level namespace
using the command

    >>> from scipy import *

Scipy sub-packages need to be imported separately, for example

    >>> from scipy import linalg, optimize


Finding Documentation
---------------------

Scipy and Numpy have HTML and PDF versions of their documentation
available at http://docs.scipy.org/, which currently details nearly
all available functionality. However, this documentation is still
work-in-progress, and some parts may be incomplete or sparse.

Python also provides the facility of documentation strings. The
functions and classes available in SciPy use this method for on-line
documentation. There are two methods for reading these messages and
getting help. Python provides the command :func:`help` in the pydoc
module. Entering this command with no arguments (i.e. ``>>> help`` )
launches an interactive help session that allows searching through the
keywords and modules available to all of Python. Running the command
help with an object as the argument displays the calling signature,
and the documentation string of the object.

The pydoc method of help is sophisticated but uses a pager to display
the text. Sometimes this can interfere with the terminal you are
running the interactive session within. A scipy-specific help system
is also available under the command scipy.info. The signature and
documentation string for the object passed to the help command are
printed to standard output (or to a writeable object passed as the
third argument). The second keyword argument of "scipy.info" defines
the maximum width of the line for printing. If a module is passed as
the argument to help than a list of the functions and classes defined
in that module is printed. For example:

.. literalinclude:: examples/1-1

Another useful command is :func:`source`. When given a function
written in Python as an argument, it prints out a listing of the
source code for that function. This can be helpful in learning about
an algorithm or understanding exactly what a function is doing with
its arguments. Also don't forget about the Python command ``dir``
which can be used to look at the namespace of a module or package.

SciPy Organization
------------------

SciPy is organized into subpackages covering different scientific
computing domains. These are summarized in the following table:

.. currentmodule:: scipy

==================  =====================================================================
Subpackage          Description
==================  =====================================================================
:mod:`cluster`      Clustering algorithms
:mod:`constants`    Physical and mathematical constants
:mod:`fftpack`      Fast Fourier Transform routines
:mod:`integrate`    Integration and ordinary differential equation solvers
:mod:`interpolate`  Interpolation and smoothing splines
:mod:`io`           Input and Output
:mod:`linalg`       Linear algebra
:mod:`maxentropy`   Maximum entropy methods
:mod:`ndimage`      N-dimensional image processing
:mod:`odr`          Orthogonal distance regression
:mod:`optimize`     Optimization and root-finding routines
:mod:`signal`       Signal processing
:mod:`sparse`       Sparse matrices and associated routines
:mod:`spatial`      Spatial data structures and algorithms
:mod:`special`      Special functions
:mod:`stats`        Statistical distributions and functions
:mod:`weave`        C/C++ integration
==================  =====================================================================

Because of their ubiquitousness, some of the functions in these
subpackages are also made available in the scipy namespace to ease
their use in interactive sessions and programs. In addition, many
basic array functions from :mod:`numpy` are also available at the
top-level of the :mod:`scipy` package. Before looking at the
sub-packages individually, we will first look at some of these common
functions.
