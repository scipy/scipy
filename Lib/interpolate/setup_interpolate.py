import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)
    
    config = default_config_dict()
    config['packages'].append(parent_package+'interpolate')
    fitpack = glob(os.path.join(local_path,'fitpack','*.f'))
    config['fortran_libraries'].append(('fitpack',{'sources':fitpack}))
    
    sources = ['_fitpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'interpolate._fitpack',sources)
    config['ext_modules'].append(ext)

    return config