============================
Low-level callback functions
============================

.. currentmodule:: scipy

Some functions in SciPy take as arguments callback functions, which
can either be python callables or low-level compiled functions.  Using
compiled callback functions can improve performance somewhat by
avoiding wrapping data in Python objects.

Such low-level functions in SciPy are wrapped in `LowLevelCallable`
objects, which can be constructed from function pointers obtained from
ctypes, cffi, Cython, or contained in Python `PyCapsule` objects.

.. autosummary::
   :toctree: generated/

   LowLevelCallable

.. seealso::

   Functions accepting low-level callables:

   `scipy.integrate.quad`, `scipy.ndimage.generic_filter`, `scipy.ndimage.generic_filter1d`,
   `scipy.ndimage.geometric_transform`

   Usage examples:

   :ref:`ndimage-ccallbacks`, :ref:`quad-callbacks`
