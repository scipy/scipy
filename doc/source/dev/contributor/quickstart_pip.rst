:orphan:

.. _quickstart-pip:

=======================================================================
Development environment quickstart guide using ``pip`` on Ubuntu Linux)
=======================================================================

This is a high-level overview of what is needed to set up a development
environment. This is only one possible way out of many. This guide assumes
you have a fresh install of Ubuntu Linux 20.04, which only has a ``python3``
executable. We also assume you have already installed ``git`` and cloned 
the SciPy repository.

Start with installing ``pip``::

    sudo apt install -y python3-pip

All further work should proceed in a virtual environment. Popular options include
the standard library ``venv`` module or a separate 
``virtualenv`` package. There are muliple third-party tutorials on how to
set up a virtual environment, so we do not cover it here. 
The rest of this guide assumes you use the ``virtualenv`` package and its
companion ``virtualenvwrapper``. (You will also need to add several
incantations to your ``.bashrc`` file.)::

    pip3 install virtualenvwrapper --user

.. note::

    We repeat: all work should happen in a virtual environment. Never use ``sudo pip``. 


Installing dependencies and building SciPy
------------------------------------------

First, you will also need the compilers for C, C++ and Fortran:: 

    sudo apt install -y gcc g++ gfortran
    
SciPy also requires BLAS and LAPACK libraries. You can install several variants
(ATLAS, OpenBLAS etc), but here we take the simplest option::

    sudo apt install -y liblapack-dev

Now, turn to the python-level work. Create a new virtual environment and activate it::

    mkvirtualenv scipy-dev

Your command prompt now lists the name of your new environment, like so
``(scipy-dev)$``. This means that the environment is active. If it is not, 
activate it manually with::

    workon scipy-dev

Inside the ``scipy-dev`` environment, install the python-level dependencies::

    pip install numpy pytest cython pybind11

Note that when the virtual environment is active, the system-wide names ``pip3``
and ``python3`` are aliased to ``pip`` and ``python``, respectively.

Now that you have all needed dependencies, build SciPy (this takes a while)::

    python setup.py build
    
And test the it::

    python runtests.py


