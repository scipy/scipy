import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict
from scipy_distutils.atlas_info import get_atlas_info

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)    
    config = default_config_dict()
    
    config['packages'].append(parent_package+'integrate')
    #config['packages'].append(parent_package+'integrate.tests')
    
    quadpack = glob(os.path.join(local_path,'quadpack','*.f'))
    config['fortran_libraries'].append(('quadpack',{'sources':quadpack}))
    
    odepack = glob(os.path.join(local_path,'odepack','*.f'))
    config['fortran_libraries'].append(('odepack',{'sources':odepack}))
    
    # should we try to weed through files and replace with calls to
    # LAPACK routines?
    linpack_lite = glob(os.path.join('integrate','linpack_lite','*.f'))
    config['fortran_libraries'].append(('linpack_lite',{'sources':linpack_lite}))
    
    mach = glob(os.path.join(local_path,'mach','*.f'))
    config['fortran_libraries'].append(('mach',{'sources':mach}))
   
    # Extension
    # flibraries.append(('blas',{'sources':blas}))
    #  Note that all extension modules will be linked against all c and
    #  fortran libraries.  But it is a good idea to at least comment
    #  the dependencies in the section for each subpackage.
    sources = ['_quadpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'integrate._quadpack',sources)
    config['ext_modules'].append(ext)
    
    # need info about blas -- how to get this???
    blas_libraries, lapack_libraries, atlas_library_dirs = get_atlas_info()
    sources = ['_odepackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(parent_package+'integrate._odepack',sources,
                    library_dirs=atlas_library_dirs,
                    libraries=['odepack', 'linpack_lite',] + blas_libraries)                    
    config['ext_modules'].append(ext)

    return config
