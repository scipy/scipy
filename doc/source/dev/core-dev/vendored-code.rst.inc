.. _vendored-code:

Vendored Code
=============
Various parts of the SciPy codebase are maintained elsewhere and vendored in SciPy.
Most of this vendored source code can be found in either the ``subprojects/`` or
``scipy/_external/`` directory.

However, some vendored source code can be found in other directories.
This is typically because code has been copied into SciPy and modified since,
or work has not yet been done to vendor that code in an idiomatic way (with
a vendoring script or git submodule).

Maintainers should be careful to *not* accept contributions
(especially trivial changes) into parts of SciPy where the code is actively maintained
upstream. Instead, they should direct contributors to the upstream repo.

Parts of the codebase to watch out for:

- ``subprojects/``
- ``scipy/_external/``
- DIRECT_, at ``scipy/optimize/_direct``
- SuperLU_, at ``scipy/sparse/linalg/_dsolve/SuperLU``
- trlib_, at ``scipy/optimize/_trlib``
- `Cython/Tempita`_, at ``scipy/_build_utils/tempita``
- `fast_matrix_market`_, at ``scipy/io/_fast_matrix_market``
- `numpydoc/docscrape`_, at ``scipy/_lib/_docscrape.py``
- `uarray`_, at ``scipy/_lib/_uarray``
- `scipy-mathjax`_, at ``doc/source/_static/scipy-mathjax``

.. _SuperLU: https://github.com/xiaoyeli/superlu
.. _trlib: https://github.com/felixlen/trlib
.. _DIRECT: https://github.com/stevengj/nlopt/tree/master/src/algs/direct
.. _`Cython/Tempita`: https://github.com/cython/cython/tree/master/Cython/Tempita
.. _`fast_matrix_market`: https://github.com/alugowski/fast_matrix_market
.. _`numpydoc/docscrape`: https://github.com/numpy/numpydoc/
.. _`uarray`: https://github.com/Quansight-Labs/uarray
.. _`scipy-mathjax`: https://github.com/scipy/scipy-mathjax

Please refer to https://github.com/scipy/scipy/issues/21232 for further details and
tracking of vendored code in the repository.
