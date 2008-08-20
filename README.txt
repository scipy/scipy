==============
fast_vectorize
==============

fast_vectorize lets you create U functions from python source code.

It works by translating the python function into C and then uses
scipy.weave to create a UFunc which calls the appropriate C function
in the inner 1-d loop.  This means that there are no Python calls
when the calculation is performed, making the calculation
fast (in particular when the arrays involved in the calculation
are very large).

Requirements:
  pypy


Use SVN to download the pypy source:
svn co http://codespeak.net/svn/pypy/dist pypy-dist

Make sure pypy can be imported, e.g. set your PYTHONPATH:
export PYTHONPATH=<path-to-pypy-dist>
