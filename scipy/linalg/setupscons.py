## Automatically adapted for scipy Oct 18, 2005 by

#!/usr/bin/env python

from __future__ import nested_scopes
import os
import sys
import re
from distutils.dep_util import newer_group, newer
from glob import glob
from os.path import join

#-------------------
# To skip wrapping single precision atlas/lapack/blas routines, set
# the following flag to True:
skip_single_routines = 0

# Some OS distributions (e.g. Redhat, Suse) provide a blas library that
# is built using incomplete blas sources that come with lapack tar-ball.
# In order to use such a library in scipy.linalg, the following flag
# must be set to True:
using_lapack_blas = 0

#--------------------

def needs_cblas_wrapper(info):
    """Returns true if needs c wrapper around cblas for calling from
    fortran."""
    import re
    r_accel = re.compile("Accelerate")
    r_vec = re.compile("vecLib")
    res = False
    try:
        tmpstr = info['extra_link_args']
        for i in tmpstr:
            if r_accel.search(i) or r_vec.search(i):
                res = True
    except KeyError:
        pass

    return res

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError

    from numpy.distutils.misc_util import Configuration

    from interface_gen import generate_interface

    config = Configuration('linalg',parent_package,top_path)

    # list of source files to register to distutils
    source_files = []

    config.add_sconscript('SConstruct')
    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(top_path='').todict())
