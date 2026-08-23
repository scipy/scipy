Building SciPy from source to use SciPy
---------------------------------------

.. note::

   If you are only trying to install SciPy, we recommend using binaries - see
   `Installation <https://scipy.org/install>`__ for details on that.

To build SciPy from source to use SciPy, there are two steps:
setting up system-level dependencies, and building SciPy itself.

.. include:: system_level.rst

Building SciPy itself
^^^^^^^^^^^^^^^^^^^^^

.. tab-set::

  .. tab-item:: Virtual env or system Python
    :sync: pip

    After all requisite system-level dependencies are installed,
    you can follow the
    `Python Packaging User Guide on creating a virtual environment <https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-virtual-environments>`__,
    then build and install SciPy via ``pip``:

    ::

      # To build the latest stable release:
      pip install scipy --no-binary scipy

      # To build a development version, you need a local clone of the SciPy git repository:
      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init
      pip install .

  .. tab-item:: Conda env
    :sync: conda

    If you are using a conda environment, ``pip`` is still the tool you use to
    invoke a from-source build of SciPy. It is important to always use the
    ``--no-build-isolation`` flag to the ``pip install`` command, to avoid
    building against a ``numpy`` wheel from PyPI. In order for that to work you
    must first install the remaining build dependencies into the conda
    environment::

      # Either install all SciPy dev dependencies into a fresh conda environment:
      conda env create -f environment.yml

      # Or, install only the required build dependencies (accurate at the time of writing):
      conda install python numpy cython pythran pybind11 compilers openblas meson-python pkg-config

      # To build the latest stable release:
      pip install scipy --no-build-isolation --no-binary scipy

      # To build a development version, you need a local clone of the SciPy git repository:
      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init
      pip install . --no-build-isolation
