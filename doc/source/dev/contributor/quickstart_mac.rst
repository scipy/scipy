.. _quickstart-mac:

================================================
Development Environment Quickstart Quide (macOS)
================================================

This quickstart guide will cover:

* setting up and maintaining a development environment, including installing compilers and SciPy dependencies;
* creating a personal fork of the SciPy repository on GitHub;
* using git to manage a local repository with development branches;
* performing an in-place build of SciPy; and 
* creating a virtual environment that adds this development version of SciPy to the Python path on macOS.

Its companion videos `Anaconda SciPy Dev: Part I (macOS)`_ and `Anaconda SciPy Dev: Part I (macOS)`_ show many of the steps being performed. This guide may diverge slightly from the videos over time with the goal of keeping this guide the simplest, up-to-date procedure.

.. note:: 

	This guide does not present the only way to set up a development environment; there are many valid choices of Python distribution, C/Fortran compiler, and installation options. The steps here can often be adapted for other choices, but we cannot provide documentation tailored for them.
	
	This guide assumes that you are starting without an existing Python 3 installation. If you already have Python 3, you might want to uninstall it first to avoid ambiguity over which Python version is being used at the command line. 

Building SciPy
--------------

(Consider following along with the companion video `Anaconda SciPy Dev: Part I (macOS)`_) 

#. Download, install, and test the latest release of the `Anaconda Distribution of Python`_. In addition to the latest version of Python 3, the Anaconda distribution includes dozens of the most popular Python packages for scientific computing, the Spyder integrated development environment (IDE), the ``conda`` package manager, and tools for managing virtual environments. 

#. Install Apple Developer Tools. An easy way to do this is to `open a terminal window <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_, enter the command ``xcode-select --install``, and follow the prompts. Apple Developer Tools includes `git <https://git-scm.com/>`_, the software we need to download and manage the SciPy source code.

#. Browse to the `SciPy repository on GitHub <https://github.com/scipy/scipy>`_ and `create your own fork <https://help.github.com/en/articles/fork-a-repo>`_. You'll need to create a GitHub account if you don't already have one.

#. Browse to your fork. Your fork will have a URL like `https://github.com/mdhaber/scipy <https://github.com/mdhaber/scipy>`_, except with your GitHub username in place of "mdhaber".

#. Click the big, green "Clone or download" button, and copy the ".git" URL to the clipboard. The URL will be the same as your fork's URL, except it will end in ".git".

#. Create a folder for the SciPy source code in a convenient place on your computer. `Navigate to it in the terminal <https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line>`_.

#. Enter the command ``git clone`` followed by your fork's .git URL. Note that this creates in the terminal's working directory a ``scipy`` folder containing the SciPy source code.

#. In the terminal, navigate into the ``scipy`` root directory (e.g. ``cd scipy``).

#. Install `Homebrew`_. Enter into the terminal |br| ``/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`` |br| or follow the installation instructions listed on the Homebrew website. Homebrew is a package manager for macOS that will help you download `gcc`, the software we will use to compile C, C++, and Fortran code included in SciPy.

#. Use Homebrew to install ``gcc`` by entering the command ``brew install gcc``.

#. In the terminal, update/upgrade all of SciPy's dependencies: ``conda update setuptools wheel cython numpy matplotlib pytest``

#. (Optional) Check your present working directory by entering ``pwd`` at the terminal. You should be in the root ``/scipy`` directory, not in a directory ending ``/scipy/scipy``.

#. Do an in-place build: enter ``python3 setup.py build_ext --inplace``. |br| This will compile the C, C++, and Fortran code that comes with SciPy. We installed ``python3`` with Anaconda. ``setup.py`` is a script in the root directory of SciPy, which is why you have to be in the SciPy root directory to call it. ``build_ext`` is a command defined in ``setup.py``, and ``--inplace`` is an option we'll use to ensure that the compiling happens in the SciPy directory you already have rather than some other folder on your computer.

#. Test the build: enter ``python3 runtests.py -v``. `runtests.py` is another script in the SciPy root directory. It runs a suite of tests that make sure SciPy is working as it should, and ``-v`` activates the ``--verbose`` option to show all the test output.

If the tests were successful, you now have a working development build of SciPy! You could stop here, but you would only be able to use this development build from within the SciPy root directory. This would be inconvenient, for instance, if you wrote a script that performs an ``import`` of something you changed in SciPy but wanted to save it elsewhere on your computer. Without taking additional steps to add this version of SciPy to the |PYTHONPATH|_ this script would ``import`` from the version of SciPy distributed with Anaconda rather than the development version you just built. (See `here <https://chrisyeh96.github.io/2017/08/08/definitive-guide-python-imports.html>`_ for much more information about how Python imports modules.)

Installing SciPy
--------------

(Consider following along with the companion video `Anaconda SciPy Dev: Part II (macOS)`_)

Currently we have *two* versions of SciPy: the latest release as installed by Anaconda, and the development version we just built. Ideally, we'd like to be able to switch between the two as needed. `Virtual environments <https://medium.freecodecamp.org/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c>`_ can do just that. With a few keystrokes in the terminal or even the click of an icon, we can enable or disable our development version. Let's set that up.

#. In a terminal window, enter ``conda list``. |br| This shows a list of all the Python packages that came with the Anaconda distribution of Python.

#. Enter ``conda create --name scipydev``. |br| This tells conda to creates a virtual environment named ``scipydev``. Note that `scipydev` can be replaced by any name you'd like to refer to your virtual environment.

#. You're still in the base environment. Activate your new virtual environment by entering ``conda activate scipydev``. |br| If you're working with an old version of ``conda``, you might need to type ``source activate scipydev`` instead (see `here <https://stackoverflow.com/questions/49600611/python-anaconda-should-i-use-conda-activate-or-source-activate-in-linux)>`_.

#. (Optional) Enter ``conda list`` again. Note that the new virtual environment has no modules installed. If you were to open a Python interpreter now, you wouldn't be able to import ``numpy``, ``scipy``, etc...

#. Enter ``conda install cython numpy matplotlib pytest spyder``. |br| Note that we're only installing SciPy's dependencies (and Spyder so we can use the IDE), but not SciPy itself.

#. Our goal now is to add our root SciPy directory to the ``PYTHONPATH`` environment variable whenever this virtual environment is activated. This will ensure that Python can find the SciPy code we are trying to `import`. This requires adding a few files and folders deep inside your Anaconda installation directory. I suggest watching `the video <https://youtu.be/Faz29u5xIZc?t=35>`_ for this part. To summarize, you want to create:

   - ``/anaconda3/envs/scipydev/conda/activate.d/env_vars.sh``, and

   - ``/anaconda3/envs/scipydev/conda/deactivate.d/env_vars.sh`` |br| where: |br|

   - ``anaconda3`` is the root directory of your Anaconda installation;

   - ``conda``, ``activate.d``, and ``deactivate.d`` are new folders; and

   - ``env_vars.sh`` is the name of two new plain text files |br| with the contents |br|
   
   - ``export PYTHONPATH=/scipy`` (where ``scipy`` is to be replaced with the full path of the SciPy root directory), and
   
   - ``unset PYTHONPATH``,
   
   respectively.

#. In a new terminal window, test your setup. If you activate your virtual environment (e.g. ``conda activate scipydev``) and run Python code that imports from SciPy, any changes you make to the SciPy code should be reflected when the code runs. After deactivating the virtual environment (``conda deactivate``), Python imports from the version of SciPy installed by Anaconda.


.. _Anaconda SciPy Dev\: Part I (macOS): https://youtu.be/1rPOSNd0ULI

.. _Anaconda SciPy Dev\: Part II (macOS): https://youtu.be/Faz29u5xIZc

.. _Anaconda Distribution of Python: https://www.anaconda.com/distribution/

.. _Homebrew: https://brew.sh/

.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables)

.. |br| raw:: html

    <br>
