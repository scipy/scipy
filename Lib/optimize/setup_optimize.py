import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)
    
    config = default_config_dict()
    config['packages'].append(parent_package+'optimize')
    minpack = glob(os.path.join(local_path,'minpack','*.f'))
    config['fortran_libraries'].append(('minpack',{'sources':minpack}))
    
    sources = ['_minpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'optimize._minpack',sources)
    config['ext_modules'].append(ext)

    return config
