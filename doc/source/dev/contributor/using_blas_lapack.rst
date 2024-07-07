.. _using-blas-lapack:

Using BLAS and LAPACK from within SciPy
=======================================

1. Cython
2. Fortran
3. C and C++


Concepts, design and common issues
----------------------------------

By design, we allow using any BLAS/LAPACK library that the user or packager
has installed on their system, and we try hard to detect or provide build
options to handle any commonly used build configuration of known BLAS and
LAPACK libraries in the wild.

When using a BLAS or LAPACK function from within SciPy, there are a number of
things we have to deal with:

1. ABI issues. BLAS and LAPACK have a Fortran ABI - and there are two different
   ABIs in use (gfortran, and f2c (or g77)). MKL and Accelerate use the f2c
   ABI, while OpenBLAS, ATLAS and Netlib BLAS *usually* use the gfortran ABI.
   However, conda-forge uses a Netlib BLAS interface with f2c ABI. Since we
   cannot know for sure which ABI we need at build time, we choose sensible
   defaults but also offer an explicit build option to override that default
   (``-Duse-g77-abi``).
2. 32-bit vs. 64-bit integers. Most libraries use 32-bit integers - this is
   called the LP64 interface. Sometimes a library is built with 64-bit integers
   (ILP64). It's also possible that both 32-bit and 64-bit symbols are present
   in the same library (this is currently the case for MKL and Accelerate). To
   accommodate ILP64, we need to use different function prototypes, changing
   the integer parameters from ``int`` to ``npy_int64``.
3. Symbol names. In general, the plain BLAS/LAPACK routine names are not enough
   to know what symbol names are present in an installed library. There may be
   a common symbol prefix and/or suffix. Examples: 

   - ``scipy-openblas`` always uses a ``scipy_`` symbol prefix.
   - The ILP64 OpenBLAS builds use a ``64_`` symbol suffix (typically;it's
     also possible that there is no suffix).
   - Accelerate uses a ``$NEWLAPACK`` suffix for LP64 symbols, and ``$NEWLAPACK_ILP64`` for ILP64 symbols.

   Finally, there's is compiler mangling. Most commonly used Fortran compilers
   (at least gfortran, Intel Fortran, and Flang) append an underscore to symbol
   names, while C compilers do not. Hence, when calling a BLAS library compiled
   with a Fortran compiler from C code, we need to append an underscore
   ourselves.

Two very useful references to learn more about SciPy's support for BLAS and
LAPACK are:

- `"Circumventing The Linker: Using SciPy's BLAS and LAPACK Within Cython" <https://proceedings.scipy.org/articles/Majora-7b98e3ed-008>`__ by Ian Henriksen (2015)
- https://pypackaging-native.github.io/key-issues/native-dependencies/blas_openmp

The most useful places to start exploring the code base are:

- ``scipy/_build_utils/src/npy_cblas.h`` for the C header that contains all the
  name mangling machinery.
- ``scipy/linalg/_generate_pyx.py`` for the machinery that generates the
  ``cython_blas``/``cython_lapack`` Cython APIs.


From Cython code
----------------

Using BLAS or LAPACK functionality from Cython code within SciPy is the most straightforward,
at least as long as you are fine with using only 32-bit (LP64) routines. If so,
then you can use ``scipy.linalg.cython_blas``/``scipy.linalg.cython_lapack`` by simply
``cimport``-ing it in your Cython code, and using the regular BLAS and LAPACK routine
names (no need to worry about symbol prefix/suffix or ABI issues).

In case you want to support 64-bit (ILP64) routines: then you need to use the
same approach as described for C/C++ below.


From Fortran code
-----------------

*TODO: describe this in a bit more detail. Not as relevant anymore perhaps,
since no one should be contributing new Fortran code anymore that needs to call
BLAS or LAPACK.*

Brief summary:

- Fortran code calls BLAS and LAPACK routines directly, without worrying about
  symbol names, compiler mangling, or ABI issues. E.g., from
  ``integrate/odepack/zvode.f``: ``CALL DZAXPY(...)``.
- Integer bitness has to be taken care of in the C wrappers for the Fortran code.
  E.g., in ``_odepackmodule.c``, ``F_INT`` is defined to ``int`` or
  ``npy_int64`` depending on whether ``HAVE_BLAS_ILP64`` is defined, and then
  function prototypes use that: ``id LSODA(lsoda_f_t *f, F_INT *neq, ...)``.
- If a symbol prefix or suffix is present, then Fortran actually really struggles.
  To deal with this, separate wrappers are added that do the name mapping (one
  source file per BLAS and LAPACK function, ~1,600 files in total - see
  `gh-19816 <https://github.com/scipy/scipy/pull/19816>`__ and
  ``scipy/_build_utils/_generate_blas_wrapper.py``).


From C/C++ code
---------------

In existing C/C++ code, there are several different approaches, usually
motivated by the same symbol names and ABI issues that we discussed higher up
on this page. We give two examples here, for SuperLU and for ``trlib``.

In SuperLU (part of ``sparse/linalg/_dsolve/``), symbol names used contain only
a trailing underscore, for example ``zgemv_``. This canonical name is mapped
within ``SuperLU/SRC/slu_Cnames.h`` to a host of other names depending on some
build config. The correct one is selected for SciPy in
``SuperLU/SRC/scipy_slu_config.h``.

In ``optimize/_trlib``, the C code uses names prefixed with ``TRLIB_``. These names are created
in ``optimize/_trlib/trlib_private.h`` - for example for ``TRLIB_DAXPY`` we have:

.. code:: c

    /* in `scipy/_build_utils/src/`, provides `BLAS_FUNC` and `CBLAS_INT` */
    #include "npy_cblas.h"

    void BLAS_FUNC(daxpy)(CBLAS_INT *n, double *alpha, double *x, CBLAS_INT *incx, double *y, CBLAS_INT *incy);

    static void trlib_daxpy(trlib_int_t *n, double *alpha, double *x, trlib_int_t *incx, double *y, trlib_int_t *incy)
    {
        CBLAS_INT n_ = *n, incx_ = *incx, incy_ = *incy;
        BLAS_FUNC(daxpy)(&n_, alpha, x, &incx_, y, &incy_);
    }

    /*
     * The lower-case to upper-case `TRLIB_DAXPY` is done with an ugly macro
     * which also uses timer functions to allow benchmarking - we leave that out
     * here (just don't do it like that!)
     */

When we're vendoring C or C++ code that is already using BLAS or LAPACK, we
typically see that that code has already figured out some way of naming their
functions (a la ``trlib_daxpy``), and we may have to adapt to it (or fix it
up).

That leaves the question what we should do when we're writing new C or C++
code. The key ingredients are ``BLAS_FUNC`` for name mangling and ``CBLAS_INT``
for LP64/ILP64 support, both defined in ``scipy/_build_utils/src/npy_cblas.h``.

**TODO finish describing canonical method here, and link to example! :**

- Practical example to fix up: https://github.com/scipy/scipy/pull/19970
  (interpolate refactor, using ``dlartg``)
- Note: CBLAS functions are already provided by the ``npy_cblas.h`` header as
  well, and NumPy uses those.
- We could use CBLAS functions in SciPy too, but don't do that as of now. Given
  that we shouldn't rely on LAPACKE, it doesn't quite make sense to rely on
  CBLAS either probably. Rather, we should use the Fortran symbols for both
  BLAS and LAPACK directly.
- To use a BLAS or LAPACK function:

  - use ``#include "scipy_lapack.h"`` (which includes ``npy_cblas.h"``)
  - if the routine you need isn't present yet, add it in ``scipy_lapack.h``.
    TBD: do we generate this header - it's basically
    ``build/scipy/linalg/_lapack_subroutines.h`` with ``C_INT`` instead of
    ``int``...
  - call the functions as ``BLAS_FUNC(name)``, where ``name`` is the canonical
    name for the BLAS or LAPACK function (e.g., ``daxpy``)


FAQ
---

**Why can't we use ``cython_blas`` from C/C++ code within SciPy?**

``scipy.linalg.cython_blas`` is a Python extension module, and not a regular shared library.
This means that the BLAS and LAPACK symbols aren't public, and hence cannot be
directly linked to. Cython knows how to get at these functions from within an
extension module; C/C++ do not. 

.. code:: zsh

   % # inspecting public symbols in the `cython_blas` extension module:
   % dyld_info -exports cython_blas.cpython-310-darwin.so
   cython_blas.cpython-310-darwin.so [arm64]:
       -exports:
           offset      symbol
           0x00000F3C  _PyInit_cython_blas

If you're now inclined to ask the follow-up question: "can't we just add such a
shared library?". Then the answer is: no, that is unlikely to work, since if
we'd provide a shared library with the plain BLAS/LAPACK symbol names, we'd get
symbol name clashes with an external BLAS library that may already be loaded by
NumPy or another package.


**Why not use ``ISO_C_BINDING``, available in Fortran 2003, to avoid ABI issues?**

For one because we are aiming to get rid of all Fortran code, so adding new Fortran
code rather than handling everything in C header files is counterproductive.
For another, we'd anyway have to handle symbol naming issues like prefix/suffixes
for ``scipy-openblas``, Accelerate, and ILP64. And it's also not magic -
Fortran itself also has two ABIs, and we also allow scenarios like using MKL
with gfortran or Netlib BLAS with Intel Fortran (both would have mismatching ABIs).


**I got everything to build, but at runtime I still get ``DLL load failed`` for ``_fblas`` - what gives?**

If you see the message::

    ImportError: DLL load failed while importing _fblas: The specified module could not be found.

without the error message telling you what the "specified module" is, what is
most likely happening is that you are on Windows and used ``python dev.py
test`` or similar, and that linked against a BLAS library that is not on the
DLL search path. To fix it, you may need to use ``dev.py
--with-scipy-openblas``, ``delvewheel``, or add the location of your shared
libraries to ``PATH``.

