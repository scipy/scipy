:orphan:

.. _quickstart-ubuntu:

=======================================================
Development environment quickstart guide (Ubuntu 16.04)
=======================================================

This quickstart guide will cover:

* setting up and maintaining a development environment, including installing compilers and SciPy build dependencies;
* creating a personal fork of the SciPy repository on GitHub;
* using git to manage a local repository with development branches;
* performing an in-place build of SciPy; and
* creating a virtual environment that adds this development version of SciPy to the Python path

in Ubuntu 16.04. *Users running Windows can follow these
instructions after setting up* `Windows Subsystem for Linux`_ *or an Amazon EC2
instance with Ubuntu 16.04 LTS. However, the instructions for setting up a
development environment with Docker may be more reliable.*

.. note::

	This guide does not present the only way to set up a development environment; there are many valid choices of Python distribution, C/Fortran compiler, and installation options. The steps here can often be adapted for other choices, but we cannot provide documentation tailored for them.

	This guide assumes that you are starting without an existing Python 3 installation. If you already have Python 3, you might want to uninstall it first to avoid ambiguity over which Python version is being used at the command line.

.. _quickstart-ubuntu-build:

Building SciPy
--------------

*Consider watching the companion videos* `Anaconda SciPy Dev: Part I (macOS)`_ and `Anaconda SciPy Dev: Part II (macOS)`_ *before starting. Although these videos are intended for macOS users, some of the procedures are the same or very similar to those described here.*

#. Download, install, and test the latest release of the `Anaconda Distribution of Python`_. In addition to the latest version of Python 3, the Anaconda distribution includes dozens of the most popular Python packages for scientific computing, the Spyder integrated development environment (IDE), the ``conda`` package manager, and tools for managing virtual environments.

   `In a terminal, navigate <https://help.ubuntu.com/community/UsingTheTerminal>`_ to the location in which you'd like to install Anaconda. You can download the file using the terminal command ``curl -O URL_OF_FILE``, where ``URL_OF_FILE`` is to be replaced with the URL of the Anaconda installer ``.sh`` file found at the Anaconda Distribution website. Run the installer by entering ``bash file.sh``, where ``file.sh`` is again to be replaced with the full name of the downloaded file. This starts the installation process. From there, simply follow the prompts, including the "Next Steps" at the end after the installer finishes.

#. (Optional) In a terminal window, enter ``conda list``. |br| This shows a list of all the Python packages that came with the Anaconda distribution of Python. Note the latest released version of SciPy is among them; this is not the development version you are going to build and will be able to modify.

   Ideally, we'd like to have both versions, and we'd like to be able to switch between the two as needed. `Virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_ can do just that. With a few keystrokes in the terminal or even the click of an icon, we can enable or disable our development version. Let's set that up.

   .. note::

      If ``conda`` is not a recognized command, try restarting your terminal. If it is still not recognized, please see "Should I add Anaconda to the macOS or Linux PATH?" in the `Anaconda FAQ`_.

#. Enter ``conda create --name scipydev``. |br| This tells ``conda`` to create a virtual environment named ``scipydev``. Note that ``scipydev`` can be replaced with any name you'd like to refer to your virtual environment.

#. You're still in the base environment. Activate your new virtual environment by entering ``conda activate scipydev``. |br| If you're working with an old version of ``conda``, you might need to type ``source activate scipydev`` instead (see `here <https://stackoverflow.com/questions/49600611/python-anaconda-should-i-use-conda-activate-or-source-activate-in-linux)>`__. Note that you'll need to have this virtual environment active whenever you want to work with the development version of SciPy.

#. (Optional) Enter ``conda list`` again. Note that the new virtual environment has no packages installed. If you were to open a Python interpreter now, you wouldn't be able to import ``numpy``, ``scipy``, etc...

#. Enter ``conda install cython numpy pytest spyder pybind11``. |br| Note that we're only installing SciPy's build dependencies (and Spyder so we can use the IDE), but not SciPy itself.

#. `Rename the file`_ ``anaconda3/envs/scipydev/lib/libgfortran.so`` to ``anaconda3/envs/scipydev/lib/libgfortran.so_backup``, where ``anaconda3`` is to be replaced with the full path of your Anaconda installation. *Note: This file provides an incorrect Fortran shared library; renaming it forces the system to find the right one for SciPy. Other libraries, however, might rely on this version, so we suggest that you only use this environment for SciPy development.*

#. Install git. An easy way to do this is to enter the command ``sudo apt install git`` in the terminal and follow the prompts. We'll use this software to download and manage the SciPy source code.

#. Browse to the `SciPy repository on GitHub <https://github.com/scipy/scipy>`_ and `create your own fork <https://help.github.com/en/articles/fork-a-repo>`_. You'll need to create a GitHub account if you don't already have one.

#. Browse to your fork. Your fork will have a URL like `https://github.com/mdhaber/scipy <https://github.com/mdhaber/scipy>`_, except with your GitHub username in place of "mdhaber".

#. Click the big, green "Clone or download" button, and copy the ".git" URL to the clipboard. The URL will be the same as your fork's URL, except it will end in ".git".

#. Create a folder for the SciPy source code in a convenient place on your computer. Navigate to it in the terminal.

#. Enter the command ``git clone`` followed by your fork's .git URL. Note that this creates in the terminal's working directory a ``scipy`` folder containing the SciPy source code.

#. In the terminal, navigate into the ``scipy`` root directory (e.g. ``cd scipy``).

#. (Optional) Check your present working directory by entering ``pwd`` at the terminal. You should be in the root ``..../scipy`` directory, not in a directory ending ``..../scipy/scipy``.

#. Install `Homebrew on Linux`_. Enter into the terminal |br| ``sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"`` |br| or follow the installation instructions listed on the `Homebrew on Linux`_ website. Homebrew requires additional packages to be installed, and lists the required commands in the terminal window. Copy and paste these commands into the terminal to complete the Homebrew setup. If no additional commands are given, they may be found `here <https://docs.brew.sh/Homebrew-on-Linux>`_.

   Homebrew is a package manager that will help you download ``gcc``, the software we will use to compile C, C++, and Fortran code included in SciPy.

#. Use Homebrew to install ``gcc`` by entering the command ``brew install gcc``.

#. Do an in-place build: enter ``python3 setup.py build_ext --inplace``. |br| This will compile the C, C++, and Fortran code that comes with SciPy. We installed ``python3`` with Anaconda. ``setup.py`` is a script in the root directory of SciPy, which is why you have to be in the SciPy root directory to call it. ``build_ext`` is a command defined in ``setup.py``, and ``--inplace`` is an option we'll use to ensure that the compiling happens in the SciPy directory you already have rather than the default location for Python packages. By building in-place, you avoid having to re-build SciPy before you can test changes to the Python code.

#. Test the build: enter ``python3 runtests.py -v``. ``runtests.py`` is another script in the SciPy root directory. It runs a suite of tests that make sure SciPy is working as it should, and ``-v`` activates the ``--verbose`` option to show all the test output. If the tests are successful, you now have a working development build of SciPy! You could stop here, but you would only be able to use this development build from within the SciPy root directory. This would be inconvenient, for instance, if you wrote a script that performs an ``import`` of something you changed in SciPy but wanted to save it elsewhere on your computer. Without taking additional steps to add this version of SciPy to the |PYTHONPATH|_ , this script would ``import`` from the version of SciPy distributed with Anaconda rather than the development version you just built. (See `here <https://chrisyeh96.github.io/2017/08/08/definitive-guide-python-imports.html>`__ for much more information about how Python imports modules.)

#. Enter ``conda develop .``, where ``.`` refers to the present directory. |br| This will allow us to ``import`` the development version of SciPy in Python regardless of Python's working directory.

#. In a new terminal window, test your setup. If you activate your virtual environment (e.g. ``conda activate scipydev``) and run Python code that imports from SciPy, any changes you make to the SciPy code should be reflected when the code runs. After deactivating the virtual environment (``conda deactivate``), Python imports from the version of SciPy installed by Anaconda. You can also check which version of SciPy you're using by executing in Python::

      import scipy
      print(scipy.__version__)

   If you have successfully imported a development version of SciPy, the word ``dev`` will appear in the output, e.g.::

      1.5.0.dev0+be97f1a


.. _Anaconda SciPy Dev\: Part I (macOS): https://youtu.be/1rPOSNd0ULI

.. _Anaconda SciPy Dev\: Part II (macOS): https://youtu.be/Faz29u5xIZc

.. _Anaconda Distribution of Python: https://www.anaconda.com/distribution/

.. _Rename the file: https://www.maketecheasier.com/rename-files-in-linux/

.. _Anaconda FAQ: https://docs.anaconda.com/anaconda/user-guide/faq/

.. _Homebrew on Linux: https://docs.brew.sh/Homebrew-on-Linux

.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install-win10

.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables

.. |br| raw:: html

    <br>
