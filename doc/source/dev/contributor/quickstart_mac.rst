:orphan:

.. _quickstart-mac:

=======================================================
Development environment quickstart guide (macOS)
=======================================================

This quickstart guide will cover:

* setting up and maintaining a development environment, including installing compilers and SciPy build dependencies;
* creating a personal fork of the SciPy repository on GitHub;
* using git to manage a local repository with development branches;
* performing an in-place build of SciPy; and
* creating a virtual environment that adds this development version of SciPy to the Python path

in macOS (Tested on 11.1).

.. note::

   This guide does not present the only way to set up a development environment; there are many valid choices of Python distribution, C/Fortran compiler, and installation options. The steps here can often be adapted for other choices, but we cannot provide documentation tailored for them.

   This guide assumes that you are starting without an existing Python 3 installation. If you already have Python 3, you might want to uninstall it first to avoid ambiguity over which Python version is being used at the command line.

.. _quickstart-mac-build:

Building SciPy
--------------

#. Install Apple Developer Tools. An easy way to do this is to `open a terminal window <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_, enter the command ``xcode-select --install``, and follow the prompts. Apple Developer Tools includes `git <https://git-scm.com/>`_, the software we need to download and manage the SciPy source code.

#. Download, install, and test the latest release of the `Anaconda Distribution of Python`_. In addition to the latest version of Python 3, the Anaconda Distribution includes dozens of the most popular Python packages for scientific computing, the ``conda`` package manager, and tools for managing virtual environments.

   If you're installing using the terminal, be sure to follow the "Next Steps"
   listed after the installer finishes. You might also need to restart your
   terminal window or enter ``source ~/.bash_profile`` for all the changes to take
   effect.

#. (Optional) In a terminal window, enter ``conda list``. This shows a list of all the Python packages that came with the Anaconda Distribution of Python. Note the latest released version of SciPy is among them; this is not the development version you are going to build and will be able to modify.

   Ideally, we'd like to have both versions, and we'd like to be able to switch between the two as needed. `Virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_ can do just that. With a few keystrokes in the terminal or even the click of an icon, we can enable or disable our development version. Let's set that up.

   .. note::

      If ``conda`` is not a recognized command, try restarting your terminal. If it is still not recognized, please see "Should I add Anaconda to the macOS or Linux PATH?" in the `Anaconda FAQ`_.

#. Enter ``conda create --name scipydev`` to create an empty virtual environment named ``scipydev`` (or another name that you prefer). This tells ``conda`` to create a new, empty environment for our packages. Activate the environment with ``conda activate scipydev``. Note that you'll need to have this virtual environment active whenever you want to work with the development version of SciPy.

#. Enter ``conda config --env --add channels conda-forge`` to tell Anaconda the source we want for our packages. Then enter ``conda install python=3.8 numpy pybind11 cython pythran pytest compilers sphinx pydata-sphinx-theme sphinx-panels matplotlib mypy`` to install the following packages:

   * ``numpy pybind11 cython pythran`` are four packages that SciPy depends on.

   * ``compilers`` holds compilers used to build SciPy's Fortran, C, and C++ source code.

   * ``pytest`` is needed for running the test suite.

   * ``sphinx``, ``pydata-sphinx-theme``, ``sphinx-panels`` and ``matplotlib`` are required to render the SciPy documentation.

   * ``mypy`` is a static type checker for Python. Consider using it.

   Note that we're installing SciPy's build dependencies and some other software, but not SciPy itself.

#. Browse to the `SciPy repository on GitHub <https://github.com/scipy/scipy>`_ and `create your own fork <https://help.github.com/en/articles/fork-a-repo>`_. You'll need to create a GitHub account if you don't already have one.

#. Browse to your fork. Your fork will have a URL like `https://github.com/mdhaber/scipy <https://github.com/mdhaber/scipy>`_, except with your GitHub username in place of "mdhaber".

#. Click the big, green "Clone or download" button, and copy the ".git" URL to the clipboard. The URL will be the same as your fork's URL, except it will end in ".git".

#. Create a folder for the SciPy source code in a convenient place on your computer. Navigate to it in the terminal.

#. Enter the command ``git clone`` followed by your fork's .git URL. Note that this creates in the terminal's working directory a ``scipy`` folder containing the SciPy source code.

#. In the terminal, navigate into the ``scipy`` root directory (e.g. ``cd scipy``).

#. Initialize git submodules: ``git submodule update --init``.

#. Do an in-place build: enter ``python3 setup.py build_ext --inplace``. |br| This will compile the C, C++, and Fortran code that comes with SciPy. We installed ``python3`` with Anaconda. ``setup.py`` is a script in the root directory of SciPy, which is why you have to be in the SciPy root directory to call it. ``build_ext`` is a command defined in ``setup.py``, and ``--inplace`` is an option we'll use to ensure that the compiling happens in the SciPy directory you already have rather than the default location for Python packages. By building in-place, you avoid having to re-build SciPy before you can test changes to the Python code.

#. Test the build: enter ``python3 runtests.py -v``. ``runtests.py`` is another script in the SciPy root directory. It runs a suite of tests that make sure SciPy is working as it should, and ``-v`` activates the ``--verbose`` option to show all the test output. If the tests are successful, you now have a working development build of SciPy! You could stop here, but you would only be able to use this development build when the Python working directory is the SciPy root directory.

#. Enter ``conda develop .``, where ``.`` refers to the present directory. |br| This will allow us to ``import`` the development version of SciPy in Python regardless of Python's working directory.

#. In a new terminal window, test your setup. If you activate your virtual environment (e.g. ``conda activate scipydev``) and run Python code that imports from SciPy, any changes you make to the SciPy code should be reflected when the code runs. After deactivating the virtual environment (``conda deactivate``), Python imports from the version of SciPy installed by Anaconda. You can also check which version of SciPy you're using by executing in Python::

      import scipy
      print(scipy.__version__)

   If you have successfully imported a development version of SciPy, the word ``dev`` will appear in the output, e.g.::

      1.6.0.dev0+be97f1a



.. _Anaconda Distribution of Python: https://www.anaconda.com/distribution/

.. _Anaconda FAQ: https://docs.anaconda.com/anaconda/user-guide/faq/


.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables

.. |br| raw:: html

    <br>
