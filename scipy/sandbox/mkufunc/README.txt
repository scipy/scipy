
mkufunc (make universal function) is a tool which lets you create
a C compiled version of a universal function (UFunc).

It works by translating the python function into C and then uses
scipy.weave to create a UFunc which calls the appropriate C function
in the inner 1-d loop.  This means that there are no Python calls
when the calculation is performed, making the calculation
fast (in particular when the arrays involved in the calculation
are very large).

Requirements:

  pypy

You need the pypy path in your PYTHONPATH environment:

$ export PYTHONPATH=/<...>/pypy-dist

