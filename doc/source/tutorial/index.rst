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
subpackage names, whereas the right column provides a description and links to the
corresponding chapter of this User Guide (if available):

.. list-table::
    :header-rows: 1

    * - Subpackage
      - Description (link to User Guide chapter)
    * - ``scipy.cluster``
      - Clustering algorithms
    * - ``scipy.constants``
      - Physical and mathematical constants
    * - ``scipy.differentiate``
      - Finite difference differentiation tools
    * - ``scipy.fft``
      - :doc:`./fft`
    * - ``scipy.fftpack``
      - Fast Fourier Transform routines (legacy)
    * - ``scipy.integrate``
      - :doc:`./integrate`
    * - ``scipy.interpolate``
      - :doc:`./interpolate`
    * - ``scipy.io``
      - :doc:`./io`
    * - ``scipy.linalg``
      - :doc:`./linalg`
    * - ``scipy.ndimage``
      - :doc:`./ndimage`
    * - ``scipy.optimize``
      - :doc:`./optimize`
    * - ``scipy.signal``
      - :doc:`./signal`
    * - ``scipy.sparse``
      - :doc:`./sparse`
    * - ``scipy.spatial``
      - :doc:`./spatial`
    * - ``scipy.special``
      - :doc:`./special`
    * - ``scipy.stats``
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

