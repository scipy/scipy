API - importing from Scipy
==========================

In Python the distinction between what is the public API of a library and what
are private implementation details is not always clear.  Unlike in other
languages like Java, it is possible in Python to access "private" function or
objects.  Occasionally this may be convenient, but be aware that if you do so
your code may break without warning in future releases.  Some widely understood
rules for what is and isn't public in Python are:

  - Methods / functions / classes and module attributes whose names begin 
    with a leading underscore are private.
  - If a class name begins with a leading underscore none of its members are 
    public, whether or not they begin with a leading underscore.
  - If a module name in a package begins with a leading underscore none of 
    its members are public, whether or not they begin with a leading
    underscore.
  - If a module or package defines ``__all__`` that authoritatively defines the
    public interface.
  - If a module or package doesn't define ``__all__`` then all names that don't 
    start with a leading underscore are public.

.. note:: Reading the above guidelines one could draw the conclusion that every
          private module or object starts with an underscore.  This is not the 
          case; the presence of underscores do mark something as private, but
          the absence of underscores do not mark something as public.

In Scipy there are modules whose names don't start with an underscore, but that
should be considered private.  To clarify which modules these are we define
below what the public API is for Scipy, and give some recommendations for how
to import modules/functions/objects from Scipy.

Guidelines for importing functions from Scipy
---------------------------------------------

The scipy namespace itself only contains functions imported from numpy.  These
functions still exist for backwards compatibility, but should be imported from
numpy directly.

Everything in the namespaces of scipy submodules is public.  In general, it is
recommended to import functions from submodule namespaces.  For example, the
function ``curve_fit`` (defined in scipy/optimize/minpack.py) should be
imported like this::

  from scipy import optimize
  result = optimize.curve_fit(...)

This form of importing submodules is preferred for all submodules except
``scipy.io`` (because ``io`` is also the name of a module in the Python
stdlib)::

  from scipy import interpolate
  from scipy import integrate
  import scipy.io as spio

In some cases, the public API is one level deeper.  For example the
``scipy.sparse.linalg`` module is public, and the functions it contains are not
available in the ``scipy.sparse`` namespace.  Sometimes it may result in more
easily understandable code if functions are imported from one level deeper.
For example, in the following it is immediately clear that ``lomax`` is a
distribution if the second form is chosen::

  # first form
  from scipy import stats
  stats.lomax(...)

  # second form
  from scipy.stats import distributions
  distributions.lomax(...)

In that case the second form can be chosen, **if** it is documented in the next
section that the submodule in question is public.


API definition
--------------

Every submodule listed below is public.  That means that these submodules are
unlikely to be renamed or changed in an incompatible way, and if that is
necessary a deprecation warning will be raised for one Scipy release before the
change is made.

* scipy.cluster

  - vq
  - hierarchy

* scipy.constants

* scipy.fftpack

* scipy.integrate

* scipy.interpolate

* scipy.io

  - arff
  - harwell_boeing
  - idl
  - matlab
  - netcdf
  - wavfile

* scipy.linalg

  - scipy.linalg.blas
  - scipy.linalg.cython_blas
  - scipy.linalg.lapack
  - scipy.linalg.cython_lapack
  - scipy.linalg.interpolative

* scipy.misc

* scipy.ndimage

* scipy.odr

* scipy.optimize

* scipy.signal

* scipy.sparse

  - linalg
  - csgraph

* scipy.spatial

  - distance

* scipy.special

* scipy.stats

  - distributions
  - mstats

* scipy.weave (deprecated)

