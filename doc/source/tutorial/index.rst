.. _user_guide:

****************
SciPy User Guide
****************

.. currentmodule:: scipy

.. sectionauthor:: Travis E. Oliphant

SciPy is a collection of mathematical algorithms and convenience functions built
on NumPy_ . It adds significant power to the interactive Python session by
providing the user with high-level commands and classes for manipulating and
visualizing data.

.. _NumPy: https://numpy.org

Subpackages
-----------

SciPy is organized into subpackages covering different scientific
computing domains. These are summarized in the following table:

==================  ======================================================
Subpackage          Description
==================  ======================================================
`cluster`           Clustering algorithms
`constants`         Physical and mathematical constants
`fftpack`           Fast Fourier Transform routines
`integrate`         Integration and ordinary differential equation solvers
`interpolate`       Interpolation and smoothing splines
`io`                Input and Output
`linalg`            Linear algebra
`ndimage`           N-dimensional image processing
`odr`               Orthogonal distance regression
`optimize`          Optimization and root-finding routines
`signal`            Signal processing
`sparse`            Sparse matrices and associated routines
`spatial`           Spatial data structures and algorithms
`special`           Special functions
`stats`             Statistical distributions and functions
==================  ======================================================

SciPy subpackages need to be imported separately, for example::

    >>> from scipy import linalg, optimize

Because of their ubiquitousness, some of the functions in these
subpackages are also made available in the `scipy` namespace to ease
their use in interactive sessions and programs. In addition, many
basic array functions from :mod:`numpy` are also available at the
top-level of the :mod:`scipy` package.

Below, you can find the complete user guide organized by subpackages.

.. raw:: latex

   \addtocontents{toc}{\protect\setcounter{tocdepth}{2}}

.. toctree::
   :caption: User guide
   :maxdepth: 1

   special
   integrate
   optimize
   interpolate
   fft
   signal
   linalg
   arpack
   csgraph
   spatial
   stats
   ndimage
   io

.. note::

   In this documentation, we will often assume that the main packages (NumPy,
   SciPy, and Matplotlib) have been imported as::

       >>> import numpy as np
       >>> import matplotlib as mpl
       >>> import matplotlib.pyplot as plt

   While we obviously don't require you to follow these conventions in your own
   code, they are highly recommended.

Below you can also find tutorials in Jupyter Notebook format.

.. toctree::
   :caption: Tutorials
   :maxdepth: 1
   
   ../notebooks/interp_transition_guide

.. raw:: latex

   \addtocontents{toc}{\protect\setcounter{tocdepth}{1}}
