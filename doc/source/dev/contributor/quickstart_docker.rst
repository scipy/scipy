:orphan:

.. _quickstart-docker:

=================================================
Development Environment Quickstart Guide (Docker)
=================================================

This document shows how to start developing SciPy in a Docker container.
These instructions should be considered a work in progress.

Docker
------

Docker is a program for running Linux virtual machines within a host
operating system. There are Docker hosts for several OS’s including
macOS, Linux, and Windows; please find the appropriate
installation instructions for your operating system at `docs.docker.com`_.

Setup a copy of the scipy repository on your computer
-----------------------------------------------------

Instructions listed here are intended for entry into a Terminal window,
or a Windows Command window. We will use *Terminal window* as a
collective term.

1) Browse to the SciPy repository on `GitHub`_ and create your `own
   fork`_. You might need to create a GitHub account if you don’t
   already have one.
2) Browse to your fork. Your fork will have a URL like
   https://github.com/andyfaff/scipy, except with your GitHub username
   in place of “andyfaff”.
3) Click the big, green “Clone or download” button, and copy the “.git”
   URL to the clipboard. The URL will be the same as your fork’s url,
   except it will end in “.git”.
4) Create a folder for the SciPy source code in a convenient place on
   your computer. `Navigate`_ to it in the terminal window.
5) Enter the command ``git clone`` followed by your fork’s .git URL.
   Note that this creates in the terminal’s working directory a
   ``scipy`` folder containing the SciPy source code. This assumes that
   you have a ``git`` command line client that is available on your
   PATH, you can follow these `instructions to install a git client`_.

Starting Docker
---------------

Instructions for getting started with Docker can be found `here`_.

1) Start up the Docker program.
2) Type ``docker info`` in your terminal window, you should see
   something like:

::

   docker info

   Containers: 0
    Running: 0
    Paused: 0
    Stopped: 0
   Images: 0
   Server Version: 17.12.0-ce
   Storage Driver: overlay2
   ...

3) In a terminal window change the directory (using the ``cd`` command)
   to the folder that contains the SciPy git repository.

4) Start up the Docker container with the following command:

::

   docker run -it --rm -v $PWD/scipy:/home/scipy scipy/scipy-dev /bin/bash

This command starts (``run``) an interactive (``-it``) Docker container
named ``scipy-dev`` (based on Ubuntu Bionic) from the ``andyfaff``
`Docker Hub repository`_. When the Docker container starts the the
``scipy`` directory from the current directory of the host (``$PWD``) is
made available in the container as ``/home/scipy``. The changes you make
from the container to any of the files in that directory are also
visible in the host, and vice-versa.

5) You should now be in the container, with something like

::

   root@468e1b9564e4:/#

as a prompt.

6) Navigate to the scipy source directory, which is shared by the host
   OS.

::

   cd /home/scipy

7) (Optional) Check your present working directory by entering ``pwd``
   at the terminal. You should be in the ``/home/scipy`` directory.

8) The container has both Python 3.6 and Python 3.7 available. To start
   using/building scipy we need to install some dependencies

::

   pip3.7 install numpy cython pytest pybind11

If you want to work on Python 3.6 use the ``pip3.6`` command.

9)  Do an in-place build. Enter
    ``python3.7 setup.py build_ext --inplace``. This will compile the C,
    C++, and Fortran code that comes with SciPy. ``setup.py`` is a
    script in the root directory of SciPy, which is why you have to be
    in the SciPy root directory to call it. ``build_ext`` is a command
    defined in ``setup.py``, and ``--inplace`` is an option we’ll use to
    ensure that the compiling happens in the SciPy directory you already
    have rather than some other folder on your computer. If you want to
    work on Python 3.6 use the ``python3.6`` command.

10) Test the build: Enter ``python3.7 runtests.py -v``. ``runtests.py``
    is another script in the SciPy root directory. It runs a suite of
    tests that make sure SciPy is working as it should, and -v activates
    the –verbose option to show all the test output.

11) You can make changes in a text editor/IDE in your host OS, and they
    will be visible within the container. Alternatively use the ``vi``
    text editor within the container to make changes. No changes made
    within the container are retained when the container is exited, only
    changes made to files/folders within mounted volumes are kept.

If you would like to contribute changes to the SciPy project, please see
:ref:`development-workflow`.

.. _NumPy: https://docs.scipy.org/doc/numpy/dev/gitwash/
.. _here: https://docs.docker.com/get-started/
.. _Docker Hub repository: https://cloud.docker.com/repository/docker/scipy/scipy-dev
.. _GitHub: https://github.com/scipy/scipy
.. _own fork: https://help.github.com/en/articles/fork-a-repo
.. _Navigate: https://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line
.. _instructions to install a git client: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _docs.docker.com: https://docs.docker.com/install/

.. |PYTHONPATH| replace:: ``PYTHONPATH``
.. _PYTHONPATH: https://docs.python.org/3/using/cmdline.html#environment-variables

.. |br| raw:: html

    <br>
