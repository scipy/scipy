============
Introduction
============

.. contents::

SciPy is a collection of mathematical algorithms and convenience
functions built on the Numpy extension of Python. It adds
significant power to the interactive Python session by providing the
user with high-level commands and classes for manipulating and
visualizing data. With SciPy an interactive Python session
becomes a data-processing and system-prototyping environment rivaling
sytems such as MATLAB, IDL, Octave, R-Lab, and SciLab.

The additional benefit of basing SciPy on Python is that this also makes a
powerful programming language available for use in developing
sophisticated programs and specialized applications. Scientific
applications using SciPy benefit from the development of
additional modules in numerous niche's of the software landscape by
developers across the world. Everything from parallel programming to
web and data-base subroutines and classes have been made available to
the Python programmer. All of this power is available in addition to
the mathematical libraries in SciPy.

This tutorial will acquaint the first-time user of SciPy with some of its most
important features. It assumes that the user has already installed the SciPy
package. Some general Python facility is also assumed, such as could be
acquired by working through the Python distribution's Tutorial. For further
introductory help the user is directed to the Numpy documentation.

For brevity and convenience, we will often assume that the main
packages (numpy, scipy, and matplotlib) have been imported as::

    >>> import numpy as np
    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt

These are the import conventions that our community has adopted
after discussion on public mailing lists.  You will see these
conventions used throughout NumPy and SciPy source code and
documentation.  While we obviously don't require you to follow
these conventions in your own code, it is highly recommended.

SciPy Organization
------------------

SciPy is organized into subpackages covering different scientific
computing domains. These are summarized in the following table:

.. currentmodule:: scipy

==================  ======================================================
Subpackage          Description
==================  ======================================================
:mod:`cluster`      Clustering algorithms
:mod:`constants`    Physical and mathematical constants
:mod:`fftpack`      Fast Fourier Transform routines
:mod:`integrate`    Integration and ordinary differential equation solvers
:mod:`interpolate`  Interpolation and smoothing splines
:mod:`io`           Input and Output
:mod:`linalg`       Linear algebra
:mod:`ndimage`      N-dimensional image processing
:mod:`odr`          Orthogonal distance regression
:mod:`optimize`     Optimization and root-finding routines
:mod:`signal`       Signal processing
:mod:`sparse`       Sparse matrices and associated routines
:mod:`spatial`      Spatial data structures and algorithms
:mod:`special`      Special functions
:mod:`stats`        Statistical distributions and functions
:mod:`weave`        C/C++ integration
==================  ======================================================

Scipy sub-packages need to be imported separately, for example::

    >>> from scipy import linalg, optimize

Because of their ubiquitousness, some of the functions in these
subpackages are also made available in the `scipy` namespace to ease
their use in interactive sessions and programs. In addition, many
basic array functions from :mod:`numpy` are also available at the
top-level of the :mod:`scipy` package. Before looking at the
sub-packages individually, we will first look at some of these common
functions.

Finding Documentation
---------------------

SciPy and NumPy have documentation versions in both HTML and PDF format
available at http://docs.scipy.org/, that cover nearly
all available functionality. However, this documentation is still
work-in-progress and some parts may be incomplete or sparse. As
we are a volunteer organization and depend on the community for
growth, your participation - everything from providing feedback to
improving the documentation and code - is welcome and actively
encouraged.

Python's documentation strings are used in SciPy for on-line
documentation. There are two methods for reading them and
getting help. One is Python's command :func:`help` in the `pydoc`
module. Entering this command with no arguments (i.e. ``>>> help`` )
launches an interactive help session that allows searching through the
keywords and modules available to all of Python. Secondly, running the command
`help(obj)` with an object as the argument displays that object's calling
signature, and documentation string.

The pydoc method of ``help`` is sophisticated but uses a pager to display
the text. Sometimes this can interfere with the terminal you are
running the interactive session within. A numpy/scipy-specific help system
is also available under the command ``numpy.info``. The signature and
documentation string for the object passed to the ``help`` command are
printed to standard output (or to a writeable object passed as the
third argument). The second keyword argument of ``numpy.info`` defines
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
