=================================================
Developing SciPy
=================================================

.. Contents::


What is SciPy?
--------------

SciPy (pronounced "Sigh Pie") is open-source software for mathematics,
science, and engineering.  It includes modules for statistics, optimization,
integration, linear algebra, Fourier transforms, signal and image processing,
ODE solvers, and more.  It is also the name of a very popular conference on
scientific programming with Python.

The SciPy library depends on NumPy, which provides convenient and fast
N-dimensional array manipulation. The SciPy library is built to work with
NumPy arrays, and provides many user-friendly and efficient numerical routines
such as routines for numerical integration and optimization. Together, they
run on all popular operating systems, are quick to install, and are free of
charge. NumPy and SciPy are easy to use, but powerful enough to be depended
upon by some of the world's leading scientists and engineers. If you need to
manipulate numbers on a computer and display or publish the results, give
SciPy a try!


SciPy structure
---------------

SciPy aims at being a robust and efficient "super-package" of a number
of modules, each of a non-trivial size and complexity.  In order for
"SciPy integration" to work flawlessly, all SciPy modules must follow
certain rules that are described in this document. Hopefully this
document will be helpful for SciPy contributors and developers as a
basic reference about the structure of the SciPy package.

Currently SciPy consists of the following files and directories:

  INSTALL.txt
    SciPy prerequisites, installation, testing, and troubleshooting.

  THANKS.txt
    SciPy developers and contributors. Please keep it up to date!!

  README.txt
    SciPy structure (this document).

  setup.py
    Script for building and installing SciPy.

  MANIFEST.in
    Additions to distutils-generated SciPy tar-balls.  Its usage is
    deprecated.

  scipy/
    Contains SciPy __init__.py and the directories of SciPy modules.

SciPy modules
+++++++++++++

In the following, a *SciPy module* is defined as a Python package, say
xxx, that is located in the scipy/ directory.  All SciPy modules should
follow the following conventions:

* Ideally, each SciPy module should be as self-contained as possible.
  That is, it should have minimal dependencies on other packages or
  modules.  Even dependencies on other SciPy modules should be kept to a
  minimum.  A dependency on NumPy is of course assumed.

* Directory ``xxx/`` must contain 

  + a file ``setup.py`` that defines
    ``configuration(parent_package='',top_path=None)`` function.  
    See below for more details.

  + a file ``info.py``. See below more details.

* Directory ``xxx/`` may contain 

  + a directory ``tests/`` that contains files ``test_<name>.py``
    corresponding to modules ``xxx/<name>{.py,.so,/}``.  See below for
    more details.

  + a file ``MANIFEST.in`` that may contain only ``include setup.py`` line.
    DO NOT specify sources in MANIFEST.in, you must specify all sources
    in setup.py file. Otherwise released SciPy tarballs will miss these sources.

  + a directory ``docs/`` for documentation.

For details, read:

  http://projects.scipy.org/numpy/wiki/DistutilsDoc


Documentation
-------------

The documentation site is here
    http://docs.scipy.org

Web sites
---------

The user's site is here
    http://www.scipy.org/

The developer's site is here
    http://projects.scipy.org/scipy/wiki


Mailing Lists
-------------

Please see the developer's list here
    http://projects.scipy.org/mailman/listinfo/scipy-dev


Bug reports
-----------

To search for bugs, please use the NIPY Bug Tracker at
    http://projects.scipy.org/scipy/query

To report a bug, please use the NIPY Bug Tracker at
    http://projects.scipy.org/scipy/newticket


License information
-------------------

See the file "LICENSE" for information on the history of this
software, terms & conditions for usage, and a DISCLAIMER OF ALL
WARRANTIES.

