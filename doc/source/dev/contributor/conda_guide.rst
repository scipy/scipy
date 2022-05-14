:orphan:

.. _conda-guide:

=====================================
Development environment guide (conda)
=====================================

The instructions below assume that you have downloaded, installed, and tested
the latest release of conda from either
`Miniforge or Mambaforge <https://github.com/conda-forge/miniforge>`_,
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
`Anaconda <https://www.anaconda.com/products/individual>`_.

   If you're installing using the terminal, be sure to follow the "Next Steps"
   listed after the installer finishes. You might also need to restart your
   terminal window or enter ``source ~/.bash_profile`` for all the changes to
   take effect.

#. (Optional) In a terminal window, enter ``conda list``. This shows a list of
   all the Python packages that came with the distribution of Python currently
   installed. Note the latest released version of SciPy may be among them; this
   is not the development version you are going to build and will be able to
   modify.

   Ideally, we'd like to have both versions, and we'd like to be able to switch
   between the two as needed. `Virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_
   can do just that. With a few keystrokes in the terminal or even the click of
   an icon, we can enable or disable our development version. Let's set that up.

   .. note::

      If ``conda`` is not a recognized command, try restarting your terminal. If
	  it is still not recognized, please see "Should I add Anaconda to the macOS
	  or Linux PATH?" in the `Anaconda FAQ`_.

#. The Python-level dependencies for building SciPy will be installed as part of
   the conda environment creation - see
   `environment.yml <https://github.com/scipy/scipy/blob/main/environment.yml>`_

   Note that we're installing SciPy's build dependencies and some other
   software, but not (yet) SciPy itself. Also note that you'll need to have
   this virtual environment active whenever you want to work with the
   development version of SciPy.

   To create the environment with all dependencies and compilers, from the root
   of the SciPy folder, do

   ::

      conda env create -f environment.yml

   Finally, activate this new environment by doing

   ::

      conda activate scipy-dev

#. Initialize git submodules: ``git submodule update --init``.

#. Do an in-place build: enter ``python setup.py build_ext --inplace``.
   This will compile the C, C++, and Fortran code that comes with SciPy. We
   installed ``python`` with Anaconda. ``setup.py`` is a script in the root
   directory of SciPy, which is why you have to be in the SciPy root directory
   to call it. ``build_ext`` is a command defined in ``setup.py``, and
   ``--inplace`` is an option we'll use to ensure that the compiling happens in
   the SciPy directory you already have rather than the default location for
   Python packages. By building in-place, you avoid having to re-build SciPy
   before you can test changes to the Python code.

#. Test the build: enter ``python runtests.py -v``.

   * ``runtests.py`` is another script in the SciPy root directory. It runs a
     suite of tests that make sure SciPy is working as it should, and ``-v``
     activates the ``--verbose`` option to show all the test output. If the
     tests are successful, you now have a working development build of SciPy!
     You could stop here, but you would only be able to use this development
     build when the Python working directory is the SciPy root directory.

#. Enter ``conda develop .``, where ``.`` refers to the present directory.
   This will allow us to ``import`` the development version of SciPy in Python
   regardless of Python's working directory.

#. In a new terminal window, or after navigating to a different folder, test
   your setup. If you activate your virtual environment (e.g.
   ``conda activate scipy-dev``) and run Python code that imports from SciPy,
   any changes you make to the SciPy code should be reflected when the code
   runs. After deactivating the virtual environment (``conda deactivate``),
   Python imports from the version of SciPy installed by your system, if any.
   You can also check which version of SciPy you're using by executing in
   Python::

      import scipy
      print(scipy.__version__)

   If you have successfully imported a development version of SciPy, the word
   ``dev`` will appear in the output, e.g.::

      1.6.0.dev0+be97f1a


.. _Anaconda Distribution of Python: https://www.anaconda.com/distribution/
.. _Anaconda FAQ: https://docs.anaconda.com/anaconda/user-guide/faq/
.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables
