pyHiGHS
=======

Cython wrappers over `HiGHS <https://github.com/ERGO-Code/HiGHS>`_.  This project is a fork of the main repo: `fork <https://github.com/mckib2/HiGHS/>`_.

The package is available on pypi as `scikit-highs <https://pypi.org/project/scikit-highs/>`_ and is tested on the following:

- Ubuntu 18, gcc 7.5.0, Python 3.6.9 and 3.7
- Windows 10, VS 2019, Python 3.8.2 x64

It has dependencies on numpy and Cython that require these to be installed before installation of scikit-highs.  Once numpy and Cython are installed, install using pip like this:

.. code:: bash

    pip install scikit-highs

Example usage from script:

.. code:: python

    from _highs.highs_wrapper import highs_wrapper
    from _highs.highs_wrapper import CONST_INF, CONST_I_INF

    # Call like this (usually using numpy arrays and A is a scipy.sparse.csc_matrix):
    res = highs_wrapper(c, A.indptr, A.indices, A.data, lhs, rhs, lb, ub, options)

`options` is a dict that supports all options of the HiGHS C++ API.

It will take a second to compile when installing, I did not create any wheels, etc, so always compiles from source.
