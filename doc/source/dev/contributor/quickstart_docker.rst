:orphan:

.. _quickstart-docker:

=================================================
Development environment quickstart guide (Docker)
=================================================

This document describes how to use a Docker container for your SciPy 
development environment.
These instructions should be considered a work in progress.

Docker
------

Docker is a program for running Linux virtual machines within a host
operating system. According to the `Docker website`_:

 A Docker container image is a lightweight, standalone, executable package of
 software that includes everything needed to run an application: code, runtime,
 system tools, system libraries and settings.
 Container images become containers at runtime, and in the case of Docker
 containers - images become containers when they run on Docker Engine.
 Available for both Linux and Windows-based applications, containerized
 software will always run the same, regardless of the host infrastructure.

Docker makes setting up a development environment easy and reliable: we
provide a Docker image with suitable compilers and Scipy's build-time 
dependencies. 
You can then use the Docker engine to execute the image as a container,  
add the latest development version of SciPy, and build SciPy.

There are Docker hosts for several OS's including:
macOS, Linux, and Windows. Please follow the appropriate
installation instructions for your operating system at `docs.docker.com`_.

.. note::

   If you have a version of an operating system that doesn't meet the
   requirements of Docker Desktop, such as Windows 10 Home,
   try `Docker Toolbox`_ .

Cloning SciPy
-------------

Before starting SciPy's Docker container, you should create a copy of the
SciPy source code on your computer. That way, you'll be able to access the
same files both from your native operating system and within the container.

.. note::
   
   Below we will use *terminal window* as a
   collective term that includes the Windows Command Prompt.

#. Browse to the `SciPy repository on GitHub`_ and `create your own fork`_.
   You'll need to create a GitHub account if you don’t
   already have one.

#. Browse to your fork. Your fork will have a URL like
   https://github.com/andyfaff/scipy, except with your GitHub username
   in place of "andyfaff".

#. Click the big, green "Clone or download" button, and copy the ".git"
   URL to the clipboard. The URL will be the same as your fork’s URL,
   except it will end in ".git".

#. Create a folder for the SciPy source code in a convenient place on
   your computer. `Navigate`_ to it in the terminal window.

#. Enter the command ``git clone`` followed by your fork’s .git URL.
   Note that this creates in the terminal’s working directory a
   ``scipy`` folder containing the SciPy source code. This assumes that
   you have a ``git`` command line client that is available on your
   PATH; if not, you can follow these `instructions to install a git client`_.

Starting Docker
---------------

Instructions for getting started with Docker can be found `here`_. After
ensuring that Docker is working correctly, follow the instructions below to
start a Docker container for SciPy development. You'll follow the same
instructions each time you want to start the container, as changes made to a
container do not persist after you close it.

#. In a terminal window, change the directory (using the ``cd`` command)
   to the root folder of the SciPy git repository, which contains the file
   ``setup.py``.

#. Ensure that Docker Desktop (or Docker Toolbox) is running, and start up the
   SciPy Docker container by entering the following command in a terminal
   window::

      docker run -it --rm -v $PWD/:/home/scipy scipy/scipy-dev:<image-tag> 
   
   If you are using Windows cmd, you may run the following command instead::

      docker run -it --rm -v %cd%:/home/scipy scipy/scipy-dev:<image-tag> 

   This command starts (``run``) an interactive (``-it``) Docker container
   named ``scipy-dev`` (based on Ubuntu focal) from the ``scipy``
   `Docker Hub repository`_. When the Docker container starts, the
   ``scipy`` directory from the current directory of the host (``$PWD``) is
   made available in the container as ``/home/scipy``. The changes you make
   from the container to any of the files in that directory are also
   visible in the host, and vice versa.

#. You should now be in the container, with something like::

      (base) root@468e1b9564e4:/home/scipy#

   as a prompt. Notice the ``(base)`` at the beginning, since we are using conda.

#. Activate the ``scipy-dev`` conda environment::

      conda activate scipy-dev

   this environment has all the dependencies you'll need to start using/building SciPy.

#. Initialize git submodules: ``git submodule update --init``.

#. Do an in-place build by entering::

      python setup.py build_ext --inplace

   This will compile the C,
   C++, and Fortran code that comes with SciPy. ``setup.py`` is a
   script in the root directory of SciPy, which is why you have to be
   in the SciPy root directory to call it. ``build_ext`` is a command
   defined in ``setup.py``, and ``--inplace`` is an option we’ll use to
   ensure that the compiling happens in the SciPy directory you already
   have rather than some other folder on your computer. 

#. Test the build by entering::

      python runtests.py -v

   ``runtests.py`` is another script in the SciPy root directory. It runs a
   suite of tests that make sure SciPy is working as it should, and ``-v``
   activates the ``–verbose`` option to show all the test output.

#. If you want to :ref:`build the documentation <rendering-documentation>`
   or import SciPy from any directory other than the SciPy root, you should
   set up SciPy for development::

      conda develop .

   where ``.`` refers to the present directory (in this case ``home/scipy``).

From here, you can start a Python console (e.g., enter ``python``) or
execute Python scripts from the command line (e.g.,
``python scriptname.py``).

You can make changes to files in the ``scipy`` directory in a text editor/IDE
in your host OS, and those changes will be reflected
within the container. Alternatively, you can use the ``vi``
text editor within the container to make changes. No changes made
within the container are retained when the container is exited; only
changes made to files/folders within mounted volumes are kept.
If you would like to contribute changes to the SciPy project, please see
:ref:`development-workflow`.

The Docker image contains all the compilers and dependencies needed for you
to work on SciPy. The Docker image uses the `Anaconda Distribution of Python`_
to manage the development environment.

.. _here: https://docs.docker.com/get-started/
.. _Docker Hub repository: https://cloud.docker.com/repository/docker/scipy/scipy-dev
.. _Scipy repository on GitHub: https://github.com/scipy/scipy
.. _create your own fork: https://help.github.com/en/articles/fork-a-repo
.. _Navigate: https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line
.. _instructions to install a git client: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _docs.docker.com: https://docs.docker.com/install/
.. _Docker website: https://www.docker.com/resources/what-container
.. _Docker Toolbox: https://docs.docker.com/toolbox/
.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables
.. _Anaconda Distribution of Python: https://www.anaconda.com/distribution/

.. |br| raw:: html

    <br>
