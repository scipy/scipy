import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)    
    config = default_config_dict()
    
    config['packages'].append(parent_package+'special')
    c_misc = glob(os.path.join(local_path,'c_misc','*.c'))
    cephes = glob(os.path.join(local_path,'cephes','*.c'))
    amos = glob(os.path.join(local_path,'amos','*.f'))
    toms = glob(os.path.join(local_path,'toms','*.f'))
    
    # C libraries
    config['libraries'].append(('c_misc',{'sources':c_misc}))
    config['libraries'].append(('cephes',{'sources':cephes}))
    
    # Fortran libraries
    config['fortran_libraries'].append(('amos',{'sources':amos}))
    config['fortran_libraries'].append(('toms',{'sources':toms}))
    
    # Extension
    # needs mach also
    sources = ['cephesmodule.c', 'amos_wrappers.c',
               'toms_wrappers.c','ufunc_extras.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'special.cephes',sources)
    config['ext_modules'].append(ext)

    return config
