.. _user_guide:

****************
SciPy User Guide
****************

.. currentmodule:: scipy

.. sectionauthor:: Travis E. Oliphant

SciPy is a collection of mathematical algorithms and convenience functions built
on NumPy_ . It adds significant power to Python by providing the user with
high-level commands and classes for manipulating and visualizing data.

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
`fft`               Discrete Fourier transforms
`fftpack`           Fast Fourier Transform routines (legacy)
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

For guidance on organizing and importing functions from SciPy subpackages, refer to the `Guidelines for Importing Functions from SciPy <https://scipy.github.io/devdocs/reference/index.html#guidelines-for-importing-functions-from-scipy>`_.

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
   sparse
   arpack
   csgraph
   spatial
   stats
   ndimage
   io


.. _executable-tutorials:

Executable tutorials
--------------------

Below you can also find tutorials in
`MyST Markdown <https://jupyterbook.org/en/stable/content/myst.html>`_ format.
These can be opened as Jupyter Notebooks with the help of the
`Jupytext <https://jupytext.readthedocs.io/en/latest/index.html>`_ extension.

.. toctree::
   :caption: Executable tutorials
   :maxdepth: 1
   
   ../notebooks/interp_transition_guide

.. raw:: latex

   \addtocontents{toc}{\protect\setcounter{tocdepth}{1}}
