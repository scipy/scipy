============
Introduction
============

.. contents::

SciPy is a collection of mathematical algorithms and convenience
functions built on the Numpy extension for Python. It adds
significant power to the interactive Python session by exposing the
user to high-level commands and classes for the manipulation and
visualization of data. With SciPy, an interactive Python session
becomes a data-processing and system-prototyping environment rivaling
sytems such as MATLAB, IDL, Octave, R-Lab, and SciLab.

The additional power of using SciPy within Python, however, is that a
powerful programming language is also available for use in developing
sophisticated programs and specialized applications. Scientific
applications written in SciPy benefit from the development of
additional modules in numerous niche's of the software landscape by
developers across the world. Everything from parallel programming to
web and data-base subroutines and classes have been made available to
the Python programmer. All of this power is available in addition to
the mathematical libraries in SciPy.

This document provides a tutorial for the first-time user of SciPy to
help get started with some of the features available in this powerful
package. It is assumed that the user has already installed the
package. Some general Python facility is also assumed such as could be
acquired by working through the Tutorial in the Python distribution.
For further introductory help the user is directed to the Numpy
documentation.

For brevity and convenience, we will often assume that the main
packages (numpy, scipy, and matplotlib) have been imported as::

    >>> import numpy as np
    >>> import scipy as sp
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
==================  ======================================================

Scipy sub-packages need to be imported separately, for example::

    >>> from scipy import linalg, optimize

Because of their ubiquitousness, some of the functions in these
subpackages are also made available in the scipy namespace to ease
their use in interactive sessions and programs. In addition, many
basic array functions from :mod:`numpy` are also available at the
top-level of the :mod:`scipy` package. Before looking at the
sub-packages individually, we will first look at some of these common
functions.

Finding Documentation
---------------------

Scipy and Numpy have HTML and PDF versions of their documentation
available at http://docs.scipy.org/, which currently details nearly
all available functionality. However, this documentation is still
work-in-progress, and some parts may be incomplete or sparse. As
we are a volunteer organization and depend on the community for
growth, your participation - everything from providing feedback to
improving the documentation and code - is welcome and actively
encouraged.

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
is also available under the command ``sp.info``. The signature and
documentation string for the object passed to the ``help`` command are
printed to standard output (or to a writeable object passed as the
third argument). The second keyword argument of ``sp.info`` defines
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
