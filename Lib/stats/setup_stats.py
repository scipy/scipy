#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package=''):
    #if parent_package:
    #    parent_package += '.'
    local_path = get_path(__name__)
    test_path = os.path.join(local_path,'tests')

    config = default_config_dict()
    #config['packages'].append(dot_join(parent_package,'stats'))
    config['packages'].append(dot_join(parent_package,'stats.tests')) 
    config['package_dir']['stats.tests'] = test_path

    # Extension
    sources = ['randmodule.c','ranlib_all.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,'stats.rand'),sources)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
