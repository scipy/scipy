#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path,default_config_dict,dot_join

def configuration(parent_package=''):
    package = 'fftpack'
    local_path = get_path(__name__)    
    config = default_config_dict(package,parent_package)

    fftpack = glob(os.path.join(local_path,'fftpack','*.f'))
    dfftpack = glob(os.path.join(local_path,'dfftpack','*.f'))

    config['fortran_libraries'].append(('fftpack',{'sources':fftpack}))
    config['fortran_libraries'].append(('dfftpack',{'sources':dfftpack}))
    
    sources = ['fftpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'fftpack'),
                    sources,
                    libraries = ['fftpack', 'dfftpack'])
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
