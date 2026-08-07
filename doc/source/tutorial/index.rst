.. _user_guide:

****************
SciPy User Guide
****************

.. currentmodule:: scipy

.. sectionauthor:: Travis E. Oliphant

SciPy is a collection of mathematical algorithms and convenience functions built
on NumPy_ . It adds significant power to Python by providing the user with
high-level commands and classes for manipulating and visualizing data.
The purpose of this user guide is to provide an overview of the functionality of
each subpackage along with some general application notes.

.. _NumPy: https://numpy.org

The following table lists the subpackages SciPy provides. The left column contains the
names which link to their :doc:`../reference/index` whereas the right column provides a
description and links to the corresponding chapter of this User Guide (if available):

.. list-table::
    :header-rows: 1

    * - Subpackage (link to API reference)
      - Description (link to User Guide chapter)
    * - :doc:`scipy.cluster <../reference/cluster>`
      - Clustering algorithms
    * - :doc:`scipy,constants <../reference/constants>`
      - Physical and mathematical constants
    * - :doc:`scipy.differentiate <../reference/differentiate>`
      - Finite difference differentiation tools
    * - :doc:`scipy.fft <../reference/fft>`
      - :doc:`./fft`
    * - :doc:`scipy.fftpack <../reference/fftpack>`
      - Fast Fourier Transform routines (legacy)
    * - :doc:`scipy.integrate <../reference/integrate>`
      - :doc:`./integrate`
    * - :doc:`scipy.interpolate <../reference/interpolate>`
      - :doc:`./interpolate`
    * - :doc:`scipy.io <../reference/io>`
      - :doc:`./io`
    * - :doc:`scipy.linalg <../reference/linalg>`
      - :doc:`./linalg`
    * - :doc:`scipy.ndimage <../reference/ndimage>`
      - :doc:`./ndimage`
    * - :doc:`scipy.optimize <../reference/optimize>`
      - :doc:`./optimize`
    * - :doc:`scipy.signal <../reference/signal>`
      - :doc:`./signal`
    * - :doc:`scipy.sparse <../reference/sparse>`
      - :doc:`./sparse`
    * - :doc:`scipy.spatial <../reference/spatial>`
      - :doc:`./spatial`
    * - :doc:`scipy.special <../reference/special>`
      - :doc:`./special`
    * - :doc:`scipy.stats <../reference/stats>`
      - :doc:`./stats`


There are also additional user guides for these topics:

- :doc:`./arpack` - Eigenvalue problem solver using iterative methods
- :doc:`./csgraph` - Compressed Sparse Graph Routines
- :ref:`scipy_parallel_execution`
- :ref:`scipy_thread_safety`
- :ref:`security`


.. toctree::
   :caption: User guide
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
   parallel_execution
   security
   thread_safety

