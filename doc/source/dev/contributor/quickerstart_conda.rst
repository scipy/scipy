.. _quickerstart-conda:

==============================================================
Development environment quickerstart guide (Linux and Mac)
==============================================================

With `Anaconda <https://www.anaconda.com/products/individual>`_ (and, on macOS,
the `Xcode Command Line Tools <https://webkit.org/build-tools/>`_) installed,
execute the following commands at the terminal from the base directory of
your `SciPy <https://github.com/scipy/scipy>`_ clone::

    # Best practice, use an environment rather than install in the base env
    conda create -n my-env
    conda activate my-env

    # Prefer conda-forge (because it has the compilers we need)
    conda config --env --add channels conda-forge

    # Now install all packages that are needed in one `conda install` command
    conda install python=3.8 numpy pybind11 cython pytest sphinx matplotlib mypy compilers

    # Build SciPy for development work plus run tests
    python runtests.py    # Alternatively, it's fine to use `python setup.py develop`

For more detailed instructions, see the other :ref:`dev-env` guides.
