#!/usr/bin/env python

import os
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path, Configuration, dot_join

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.system_info import get_info
    package = 'signal'
    local_path = get_path(__name__,parent_path)    
    config = Configuration(package, parent_package)

    config.add_extension('sigtools',
        sources=['sigtoolsmodule.c','firfilter.c','medianfilter.c'],
    )
    
    config.add_extension('spline',
        sources = ['splinemodule.c','S_bspline_util.c','D_bspline_util.c',
                   'C_bspline_util.c','Z_bspline_util.c','bspline_util.c'],
    )
    
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(**configuration(parent_path=''))
