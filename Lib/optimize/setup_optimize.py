#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package=''):
    config = default_config_dict('optimize',parent_package)
    local_path = get_path(__name__)    

    if parent_package:
        config['packages'].append(parent_package+'.optimize')
 
    test_path = os.path.join(local_path,'tests')
    config['packages'].append(dot_join(parent_package,'optimize.tests'))
    config['package_dir']['optimize.tests'] = test_path       

    minpack = glob(os.path.join(local_path,'minpack','*.f'))
    config['fortran_libraries'].append(('minpack',{'sources':minpack}))
    
    sources = ['_minpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'optimize._minpack',sources,
                    libraries = ['minpack'])
    config['ext_modules'].append(ext)

    rootfind = glob(os.path.join(local_path,'Zeros','*.c'))
    roothead = os.path.join(local_path,'zeros.h')
    config['libraries'].append(('rootfind',{'sources':rootfind,
                                            'headers':roothead}))
    sources = ['zeros.c']
    sources = [os.path.join(local_path,x) for x in sources]    
    ext = Extension(dot_join(parent_package,"optimize._zeros"),sources,
                    libraries=['rootfind'])
    config['ext_modules'].append(ext)
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
