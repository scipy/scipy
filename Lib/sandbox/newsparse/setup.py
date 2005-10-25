#!/usr/bin/env python

# EJS: running 'setup.py install' doesn't yet install files to the right place.
# The basic spmatrixmodule.c compiles okay on my machine, but I haven't yet
# tested it.


import os
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path,Configuration,dot_join
join = os.path.join
import glob

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.system_info import get_info
    config = Configuration('pysparse', parent_package, parent_path)
    #local_path = get_path(__name__,parent_path)

    config.add_data_dir('docs')
    config.add_data_dir('examples')
    config.add_data_dir('tests')
    #config.add_data_files(('include'))
    headers = glob.glob(os.path.join ("include","pysparse","*.h"))
    #header_dir = join(*(config.name.split('.')+['include','scipy']))
    #headers = ["csr_mat.h", "spmatrix.h", "spmatrix_api.h", "sss_mat.h", "ll_mat.h", "spmatrix.h"]
    #config.add_headers(headers)
    config.add_extension('pysparse',
        sources = ['src/spmatrixmodule.c'],
        include_dirs = ['include/']
        )

    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(description="pysparse - A port of Roman Geus's PySparse to newscipy",
          author="Ed Schofield, Robert Cimrman",
          author_email = "edschofield@scipy.org, cimrman3@ntc.zcu.cz",
          maintainer_email = "scipy-dev@scipy.org",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration(parent_path='').todict())

################


