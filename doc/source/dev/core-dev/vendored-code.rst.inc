.. _vendored-code:

Vendored Code
=============
Many parts of the SciPy codebase are maintained elsewhere, and vendored in SciPy.
Some of these parts are vendored as git submodules, for example, ``boost_math``.

Other parts are not vendored as git submodules, despite having a maintained upstream.
This is mainly for historical reasons, and it is possible that some of these parts
will see patches contributed upstream and become git submodules in the future.

Maintainers should be careful to *not* accept contributions
(especially trivial changes) into parts of SciPy where the code is actively maintained
upstream. Instead, they should direct contributors to the upstream repo.
Currently, this includes the following parts of the codebase:

- DIRECT_, at ``scipy/optimize/_direct``
- ARPACK_, at ``scipy/sparse/linalg/_eigen/arpack/ARPACK``
- SuperLU_, at ``scipy/sparse/linalg/_dsolve/SuperLU``
- QHull_, at ``scipy/spatial/qhull_src``
- trlib_, at ``scipy/optimize/_trlib``
- UNU.RAN_, at ``scipy/stats/_unuran``

.. _ARPACK: https://github.com/opencollab/arpack-ng
.. _SuperLU: https://github.com/xiaoyeli/superlu
.. _QHull: https://github.com/qhull/qhull
.. _trlib: https://github.com/felixlen/trlib
.. _UNU.RAN: https://statmath.wu.ac.at/unuran/
.. _DIRECT: https://github.com/stevengj/nlopt/tree/master/src/algs/direct
