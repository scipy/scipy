:orphan:

.. _quickstart-pip:

======================================================================
Development environment quickstart guide using ``pip`` on Ubuntu Linux
======================================================================

This is a high-level overview of what is needed to set up a development
environment. This is only one possible way out of many. This guide assumes
you have a fresh install of Ubuntu Linux 20.04, which only has a ``python3``
executable. We also assume you have already installed ``git`` and cloned
the SciPy repository.


Installing the system-level dependencies
----------------------------------------

First, you will also need the compilers for C, C++ and Fortran::

    sudo apt install -y gcc g++ gfortran

SciPy also requires BLAS and LAPACK libraries. You can install several variants
(ATLAS, OpenBLAS etc), but here we take the simplest option::

    sudo apt install -y liblapack-dev


Installing the python-level dependencies
----------------------------------------

Start with installing ``pip``::

    sudo apt install -y python3-pip

All further work should proceed in a virtual environment. Popular options include
the standard library ``venv`` module or a separate
``virtualenv`` package. There are muliple third-party tutorials on how to
set up a virtual environment, so we cover only briefly these two options
here.

.. note::

    We repeat: all work should happen in a virtual environment. Never use ``sudo pip``.


Using ``virtualenv``
~~~~~~~~~~~~~~~~~~~~

Install the ``virtualenvwrapper`` package::

    python3 -m pip install virtualenvwrapper --user

Edit the ``.bashrc`` file to add some environment variables which are used
internally by the ``virtualenvwrapper``::

    export WORKON_HOME=$HOME/virtualenvs
    export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3
    . $HOME/.local/bin/virtualenvwrapper.sh

Here we store the virtualenvs in a ``virtualenvs`` folder in the home directory.
(you might need to create the folder manually).

Now open a new terminal window for the changes to the ``.bashrc`` to take effect.

Create a new virtual environment and activate it::

    mkvirtualenv scipy-dev

Your command prompt now lists the name of your new environment, like so
``(scipy-dev)$``. This means that the environment is active. If it is not,
activate it manually with::

    workon scipy-dev

Note ``mkvirtualenv`` and ``workon`` commands come from the ``virtualwrapper``
package.



Using the standard-library ``venv`` package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the ``venv`` package::

    sudo apt install -y python3-venv

Change the directory to your home folder and create a directory ``.venvs`` there.
Create the virtualenvironment::

    python3 -m venv scipy-dev

To activate the environment, use ::

    source $HOME/.venvs/scipy-dev/bin/activate

Your command prompt now lists the name of your new environment, like so
``(scipy-dev)$``.

(For the official docs for the ``venv`` package see
https://docs.python.org/3/tutorial/venv.html).


Building SciPy
--------------

Inside the ``scipy-dev`` environment, install the python-level dependencies::

    python -m pip install numpy pytest cython pythran pybind11

Note that when the virtual environment is active, the system-wide names ``pip3``
and ``python3`` are aliased to ``pip`` and ``python``, respectively.

Now that you have all needed dependencies, navigate to the directory where
you cloned the source code into, and build SciPy (this takes a while)::

    python setup.py build

Optionally, test it::

    python runtests.py


