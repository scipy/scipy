#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package='',parent_path=None):
    package_name = 'interpolate'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package_name, parent_package)

    fitpack = glob(os.path.join(local_path,'fitpack','*.f'))
    config['fortran_libraries'].append(('fitpack',{'sources':fitpack}))
    
    sources = ['_fitpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package_name,'_fitpack'),
                    sources,
                    libraries = ['fitpack'])
    config['ext_modules'].append(ext)

    # New interface to fitpack, requires f2py with usercode support:
    sources = ['fitpack.pyf']
    sources = [os.path.join(local_path,x) for x in sources]
    ext_args = {
        'name': dot_join(parent_package,package_name,'dfitpack'),
        'sources': sources,
        'libraries': ['fitpack'],
        }
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
