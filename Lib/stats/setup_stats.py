#!/usr/bin/env python

import os
from glob import glob
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path, Configuration, dot_join

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.system_info import get_info, dict_append
    package = 'stats'
    local_path = get_path(__name__,parent_path)
    config = Configuration(package, parent_package)

    statlib = glob(os.path.join(local_path, 'statlib','*.f'))

    config.add_library('statlib',
                       sources=statlib)
    
    # add statlib module
    config.add_extension('statlib',
        sources=['statlib.pyf'],
        f2py_options=['--no-wrap-functions'],
        libraries=['statlib'],
    )
    
    # add futil module
    config.add_extension('futil',
        sources=['futil.f'],
    )

    # add mvn module
    config.add_extension('mvn',
        sources=['mvn.pyf','mvndst.f'],
    )

    return config

if __name__ == '__main__':    
    from scipy.distutils.core import setup
    setup(**configuration(parent_path=''))
