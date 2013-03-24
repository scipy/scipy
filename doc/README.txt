SciPy Documentation
===================

How to build it
---------------
The easy way to build the documentation is to run

    python setup.py build_sphinx

This will first build Scipy in-place, and then generate documentation for it.

Another way
-----------
1. Optionally download an XML dump of the newest docstrings from the doc wiki
   at ``/pydocweb/dump`` and save it as ``dump.xml``.
2. Run ``make html`` or ``make dist``

Note that ``make html`` builds the documentation for the currently installed
version of Scipy, not the one corresponding to the source code here.

