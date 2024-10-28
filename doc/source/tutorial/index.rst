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

Subpackages and User Guides
---------------------------

SciPy is organized into subpackages covering different scientific
computing domains. These are summarized in the following table, with
their API reference linked in the Subpackage column, and user guide (if available)
linked in the Description column:

==================  ========================================
Subpackage          Description and User Guide
==================  ========================================
`cluster`           Clustering algorithms
`constants`         Physical and mathematical constants
`differentiate`     Finite difference differentiation tools
`fft`               :doc:`./fft`
`fftpack`           Fast Fourier Transform routines (legacy)
`integrate`         :doc:`./integrate`
`interpolate`       :doc:`./interpolate`
`io`                :doc:`./io`
`linalg`            :doc:`./linalg`
`ndimage`           :doc:`./ndimage`
`odr`               Orthogonal distance regression
`optimize`          :doc:`./optimize`
`signal`            :doc:`./signal`
`sparse`            :doc:`./sparse`
`spatial`           :doc:`./spatial`
`special`           :doc:`./special`
`stats`             :doc:`./stats`
==================  ========================================

There are also additional user guides for these topics:

- :doc:`./arpack` - Eigenvalue problem solver using iterative methods
- :doc:`./csgraph` - Compressed Sparse Graph Routines

For guidance on organizing and importing functions from SciPy subpackages, refer to the `Guidelines for Importing Functions from SciPy <https://scipy.github.io/devdocs/reference/index.html#guidelines-for-importing-functions-from-scipy>`_.

.. raw:: latex

   \addtocontents{toc}{\protect\setcounter{tocdepth}{2}}

.. toctree::
   :caption: User guide
   :maxdepth: 1
   :hidden:

   fft
   integrate
   interpolate
   io
   linalg
   ndimage
   optimize
   signal
   sparse
   spatial
   special
   stats
   arpack
   csgraph

.. raw:: latex

   \addtocontents{toc}{\protect\setcounter{tocdepth}{1}}
