.. _vendored-code:

Vendored Code
=============
Many parts of the SciPy codebase are maintained elsewhere, and vendored in SciPy.
Some of these parts are vendored as git submodules, for example, ``xsf``,
or are placed under the ``subprojects`` directory (or both).

Other parts are not vendored as git submodules or under the ``subprojects`` directory,
despite having a maintained upstream. This is usually either because:

1. a subset of the upstream repo is vendored with a script.
    (It is possible that these parts will be moved under ``subprojects`` in the future.)

2. Code has been copied into SciPy and modified since.
    (It is possible that some of these parts will see patches contributed upstream and
    become git submodules or be moved under ``subprojects`` in the future.)

Maintainers should be careful to *not* accept contributions
(especially trivial changes) into parts of SciPy where the code is actively maintained
upstream. Instead, they should direct contributors to the upstream repo.

Parts of the codebase which are vendored with a script include:

- PRIMA_, at ``scipy/_lib/pyprima``

Parts of the codebase which contain code copied from an upstream include:

- DIRECT_, at ``scipy/optimize/_direct``
- ARPACK_, at ``scipy/sparse/linalg/_eigen/arpack/ARPACK``
- SuperLU_, at ``scipy/sparse/linalg/_dsolve/SuperLU``
- trlib_, at ``scipy/optimize/_trlib``
- UNU.RAN_, at ``scipy/stats/_unuran``
- `Cython/Tempita`_, at ``scipy/_build_utils/tempita``
- `fast_matrix_market`_, at ``scipy/io/_fast_matrix_market``
- `numpydoc/docscrape`_, at ``scipy/_lib/_docscrape.py``

.. _ARPACK: https://github.com/opencollab/arpack-ng
.. _SuperLU: https://github.com/xiaoyeli/superlu
.. _trlib: https://github.com/felixlen/trlib
.. _UNU.RAN: https://statmath.wu.ac.at/unuran/
.. _DIRECT: https://github.com/stevengj/nlopt/tree/master/src/algs/direct
.. _PRIMA: https://github.com/libprima/prima
.. _`Cython/Tempita`: https://github.com/cython/cython/tree/master/Cython/Tempita
.. _`fast_matrix_market`: https://github.com/alugowski/fast_matrix_market
.. _`numpydoc/docscrape`: https://github.com/numpy/numpydoc/

Please refer to https://github.com/scipy/scipy/issues/21232 for further details and
tracking of vendored code in the repository.
