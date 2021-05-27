.. _quickerstart-conda:

==============================================================
Development environment quickerstart guide (Linux and Mac)
==============================================================

With conda installed (through `Miniforge or Mambaforge <https://github.com/conda-forge/miniforge>`_,
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
`Anaconda <https://www.anaconda.com/products/individual>`_),
execute the following commands at the terminal from the base directory of
your `SciPy <https://github.com/scipy/scipy>`_ clone::

    # Create an environment with all development dependencies
    conda env create -f environment.yml  # works with `mamba` too
    conda activate scipy-dev

    # Initialize git submodules
    git submodule update --init

    # Build SciPy for development work plus run tests
    python runtests.py    # Alternatively, it's fine to use `python setup.py develop`

    # Install SciPy in develop mode in conda
    conda develop .

For more detailed instructions, see the other :ref:`dev-env` guides.
