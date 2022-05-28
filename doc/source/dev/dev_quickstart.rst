.. _dev-quickstart:

============================
Contributor quickstart guide
============================

After :ref:`getting the source code from GitHub <git-start>`, there are three
steps to start contributing:

1. **Set up a development environment**

   Using ``conda``, or some flavor of the many virtual environment management
   tools, you can make sure the development version of SciPy does not interfere
   with any other local installations of SciPy on your machine.

2. **Build SciPy**

   SciPy uses compiled code for speed, which means you might need extra
   dependencies to complete this step depending on your system.

3. **Perform development tasks**

   These can include any changes you want to make to the source code, running
   tests, building the documentation, running benchmarks, etc.

Basic workflow
==============

.. note::

    We **strongly** recommend using a virtual environment setup, such as
    ``venv`` or ``conda``.

Since SciPy contains parts written in C, C++, and Fortran that need to be
compiled before use, make sure you have the necessary compilers and Python
development headers installed. If you are using ``conda``, these will be
installed automatically. If you are using ``pip``, check which
:ref:`system-level dependencies <system-level>` you might need.

First, fork a copy of the main SciPy repository in GitHub onto your own
account and then create your local repository via::

    git clone git@github.com:YOURUSERNAME/scipy.git scipy
    cd scipy
    git submodule update --init
    git remote add upstream https://github.com/scipy/scipy.git

Next, set up your development environment.

.. tabs::

    .. tab:: conda

        With ``conda`` installed (through
        `Miniforge or Mambaforge <https://github.com/conda-forge/miniforge>`_,
        `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
        `Anaconda <https://www.anaconda.com/products/individual>`_), execute the
        following commands at the terminal from the base directory of your
        `SciPy <https://github.com/scipy/scipy>`_ clone::

            # Create an environment with all development dependencies
            conda env create -f environment.yml  # works with `mamba` too
            # Activate the environment
            conda activate scipy-dev

        Your command prompt now lists the name of your new environment, like so
        ``(scipy-dev)$``.

    .. tab:: pip+venv

        **With** :ref:`system-level dependencies <system-level>` **installed**, execute
        the following commands at the terminal from the base directory of your
        `SciPy <https://github.com/scipy/scipy>`_ clone:

        .. code:: bash

            # Create the virtual environment
            python -m venv scipy-dev
            # Activate the environment
            source $HOME/.venvs/scipy-dev/bin/activate
            # Install python-level dependencies
            python -m pip install numpy pytest cython pythran pybind11 meson ninja

        Your command prompt now lists the name of your new environment, like so
        ``(scipy-dev)$``.

Finally, build SciPy for development and run the test suite to make sure your
installation is successful. On Linux and OSX, you should use::

    python dev.py

This builds SciPy first, so the first time it may take some time.

If you run into a build issue, or need more detailed build documentation
including building on Windows, see :ref:`building`.

Some of the tests in SciPy are very slow and need to be separately
enabled. See :ref:`runtests` for details.

Other workflows
===============

This is only one possible way to set up your development environment out of
many. For more detailed instructions, see the :ref:`contributor-toc`.
