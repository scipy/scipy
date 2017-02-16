==================================
 Building an OSX binary for scipy
==================================

This directory contains the scripts to build a universal binary for
OSX.  The binaries work on OSX 10.4 and 10.5.

The docstring in build.py may contain more current details.

Requirements
============

* bdist_mpkg v0.4.3

Build
=====

The build script will build a scipy distribution using bdist_mpkg and
create the mac package (mpkg) bundled in a disk image (dmg).  To run
the build script::

  python build.py

Install and test
----------------

The *install_and_test.py* script will find the scipy*.mpkg, install it
using the Mac installer and then run the scipy test suite.  To run the
install and test::

  python install_and_test.py

