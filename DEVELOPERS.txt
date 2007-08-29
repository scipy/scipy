.. -*- rest -*-
.. NB! Keep this document a valid restructured document.


Developing SciPy
================

:Author: Pearu Peterson <pearu@cens.ioc.ee>
:Modified by: Ed Schofield <edschofield@gmail.com>
:Last changed: $Date$
:Revision: $Revision$
:Discussions to: scipy-dev@scipy.org

.. Contents::

Introduction
------------

SciPy aims at being a robust and efficient "super-package" of a number
of modules, each of a non-trivial size and complexity.  In order for
"SciPy integration" to work flawlessly, all SciPy modules must follow
certain rules that are described in this document. Hopefully this
document will be helpful for SciPy contributors and developers as a
basic reference about the structure of the SciPy package.

SciPy structure
---------------

Currently SciPy consists of the following files and directories:

  INSTALL.txt
    SciPy prerequisites, installation, testing, and troubleshooting.

  THANKS.txt
    SciPy developers and contributors. Please keep it up to date!!

  DEVELOPERS.txt
    SciPy structure (this document).

  setup.py
    Script for building and installing SciPy.

  MANIFEST.in
    Additions to distutils-generated SciPy tar-balls.  Its usage is
    deprecated.

  scipy/
    Contains SciPy __init__.py and the directories of SciPy modules.




SciPy modules
-------------

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

  http://svn.scipy.org/svn/numpy/trunk/numpy/doc/DISTUTILS.txt

Open issues and discussion
--------------------------

Documentation
+++++++++++++

This is an important feature where SciPy is currently lacking. A few
SciPy modules have some documentation but they use different formats
and are mostly out of date.  We could use some help with this.

Currently there are

* A SciPy tutorial by Travis E. Oliphant.  This is maintained using LyX. 
  The main advantage of this approach is that one can use mathematical
  formulas in documentation.

* I (Pearu) have used reStructuredText formated .txt files to document
  various bits of software. This is mainly because ``docutils`` might
  become a standard tool to document Python modules. The disadvantage
  is that it does not support mathematical formulas (though, we might
  add this feature ourself using e.g. LaTeX syntax).

* Various text files with almost no formatting and mostly badly out
  dated.

* Documentation strings of Python functions, classes, and modules.
  Some SciPy modules are well-documented in this sense, others are very
  poorly documented. Another issue is that there is no consensus on how
  to format documentation strings, mainly because we haven't decided
  which tool to use to generate, for instance, HTML pages of
  documentation strings.

So, we need unique rules for documenting SciPy modules. Here are some
requirements that documentation tools should satsify:

* Easy to use. This is important to lower the threshold of developers
  to use the same documentation utilities.

* In general, all functions that are visible to SciPy end-users, must
  have well-maintained documentation strings.

* Support for mathematical formulas. Since SciPy is a tool for
  scientific work, it is hard to avoid formulas to describe how its
  modules are good for. So, documentation tools should support LaTeX.

* Documentation of a feature should be closely related to its
  interface and implementation. This is important for keeping
  documentation up to date. One option would be to maintain
  documentation in source files (and have a tool that extracts
  documentation from sources). The main disadvantage with that is the
  lack of convenience writing documentation as the editor would be in
  different mode (e.g. Python mode) from the mode suitable for
  documentation.

* Differentiation of implementation (e.g. from scanning sources) and
  concept (e.g. tutorial, users guide, manual) based docs.
  


