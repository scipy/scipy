#!/usr/bin/env python

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#    Python distutils setup for PyNifti
#
#    Copyright (C) 2006-2007 by
#    Michael Hanke <michael.hanke@gmail.com>
#
#    This is free software; you can redistribute it and/or
#    modify it under the terms of the MIT License.
#
#    This package is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the COPYING
#    file that comes with this package for more details.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

from distutils.core import setup, Extension
import os
import numpy
from glob import glob

nifti_wrapper_file = os.path.join('nifti', 'nifticlib.py')

# create an empty file to workaround crappy swig wrapper installation
if not os.path.isfile(nifti_wrapper_file):
    open(nifti_wrapper_file, 'w')

# find numpy headers
numpy_headers = os.path.join(os.path.dirname(numpy.__file__),'core','include')


# Notes on the setup
# Version scheme is:
# 0.<4-digit-year><2-digit-month><2-digit-day>.<ever-increasing-integer>

setup(name       = 'pynifti',
    version      = '0.20070930.1',
    author       = 'Michael Hanke',
    author_email = 'michael.hanke@gmail.com',
    license      = 'MIT License',
    url          = 'http://apsy.gse.uni-magdeburg.de/hanke',
    description  = 'Python interface for the NIfTI IO libraries',
    long_description = """ """,
    packages     = [ 'nifti' ],
    scripts      = glob( 'bin/*' ),
    ext_modules  = [ Extension( 'nifti._nifticlib', [ 'nifti/nifticlib.i' ],
            include_dirs = [ '/usr/include/nifti', numpy_headers ],
            libraries    = [ 'niftiio' ],
            swig_opts    = [ '-I/usr/include/nifti',
                             '-I' + numpy_headers ] ) ]
    )

