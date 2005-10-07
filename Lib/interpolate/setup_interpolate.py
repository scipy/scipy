#!/usr/bin/env python

import os
from glob import glob
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path, Configuration, dot_join

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.system_info import get_info, dict_append
    package_name = 'interpolate'
    local_path = get_path(__name__,parent_path)
    config = Configuration(package_name, parent_package)

    config.add_library('fitpack',
                       sources=[os.path.join('fitpack', '*.f')],
                      )
                      
    config.add_extension('_fitpack',
                         sources=['_fitpackmodule.c'],
                         libraries=['fitpack'],
                        )

    config.add_extension('dfitpack',
                         sources=['fitpack.pyf'],
                         libraries=['fitpack'],
                        )

    return config

if __name__ == '__main__':    
    from scipy.distutils.core import setup
    setup(**configuration(parent_path=''))
