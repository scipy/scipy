#!/usr/bin/env python

import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package='',parent_path=None):
    package = 'signal'
    local_path = get_path(__name__,parent_path)    
    config = default_config_dict(package, parent_package)

    sources = ['sigtoolsmodule.c','firfilter.c','medianfilter.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'sigtools'), sources)
    config['ext_modules'].append(ext)
    
    sources = ['splinemodule.c','S_bspline_util.c','D_bspline_util.c',
               'C_bspline_util.c','Z_bspline_util.c','bspline_util.c']
    sources = [os.path.join(local_path,x) for x in sources]               
    ext = Extension(dot_join(parent_package,package,'spline'),sources)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
