#!/usr/bin/env python

"""Note: to be able to import this package after installation you need to
add a symbolic link manually in the site-packages/scipy/sandbox/pysparse/
directory from pysparse.so to spmatrix.so.  Then you can import using:

    from scipy.sandbox import spmatrix

(The reason for this is to keep the patch set from PySparse minimal while
working around PySparse's weird inconsistency in its module name.)
"""

from os.path import join
from glob import glob


def configuration(parent_package='',parent_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('pysparse', parent_package, parent_path)

    config.add_data_dir('docs')
    config.add_data_dir('examples')
    config.add_data_dir('tests')
    headers = glob(join("include","pysparse","*.h"))
    config.add_extension('pysparse',
        sources = ['src/spmatrixmodule.c'],
        include_dirs = ['include/']
        )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(description="pysparse - A port of Roman Geus's PySparse to numpy",
          author="Ed Schofield, Robert Cimrman",
          author_email = "edschofield@scipy.org, cimrman3@ntc.zcu.cz",
          maintainer_email = "scipy-dev@scipy.org",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration(parent_path='').todict())
