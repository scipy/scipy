#!/usr/bin/env python

import os
import sys
from glob import glob
import shutil

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import get_path,\
         default_config_dict, dot_join
    from scipy_distutils.system_info import dict_append

    package = 'special'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)

    define_macros = []
    if sys.byteorder == "little":
        define_macros.append(('USE_MCONF_LE',None))
    else:
        define_macros.append(('USE_MCONF_BE',None))
    if sys.platform=='win32':
        define_macros.append(('NOINFINITIES',None))
        define_macros.append(('NONANS',None))

    c_misc = glob(os.path.join(local_path,'c_misc','*.c'))
    cephes = glob(os.path.join(local_path,'cephes','*.c'))
    if sys.platform=='win32':
        cephes = [f for f in cephes if os.path.basename(f)!='fabs.c']
    mach = glob(os.path.join(local_path,'mach','*.f'))
    amos = glob(os.path.join(local_path,'amos','*.f'))
    toms = glob(os.path.join(local_path,'toms','*.f'))
    cdf = glob(os.path.join(local_path,'cdflib','*.f'))
    specfun = glob(os.path.join(local_path, 'specfun','*.f'))
    
    # C libraries
    config['libraries'].append(('c_misc',{'sources':c_misc}))
    config['libraries'].append(('cephes',{'sources':cephes,
                                          'macros':define_macros}))

    # Fortran libraries
    config['fortran_libraries'].append(('mach',{'sources':mach}))
    config['fortran_libraries'].append(('amos',{'sources':amos}))
    config['fortran_libraries'].append(('toms',{'sources':toms}))
    config['fortran_libraries'].append(('cdf',{'sources':cdf}))
    config['fortran_libraries'].append(('specfun',{'sources':specfun}))
    
    # Extension
    sources = ['cephesmodule.c', 'amos_wrappers.c', 'specfun_wrappers.c',
               'toms_wrappers.c','cdf_wrappers.c','ufunc_extras.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'cephes'),sources,
                    libraries = ['amos','toms','c_misc','cephes','mach',
                                 'cdf', 'specfun'],
                    define_macros = define_macros
                    )
    config['ext_modules'].append(ext)

    ext_args = {'name':dot_join(parent_package,package,'specfun'),
                'sources':[os.path.join(local_path,'specfun.pyf')],
                'f2py_options':['--no-wrap-functions'],
                #'define_macros':[('F2PY_REPORT_ATEXIT_DISABLE',None)],
                'libraries' : ['specfun']
                }
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    return config

def get_package_config(name):
    sys.path.insert(0,os.path.join('scipy_core',name))
    try:
        mod = __import__('setup_'+name)
        config = mod.configuration()
    finally:
        del sys.path[0]
    return config

if __name__ == '__main__':
    extra_packages = []
    try: import scipy_base
    except ImportError: extra_packages.append('scipy_base')
    try: import scipy_test
    except ImportError: extra_packages.append('scipy_test')
    try: import scipy_distutils
    except ImportError:
        extra_packages.append('scipy_distutils')
        sys.path.insert(0,'scipy_core')

    from scipy_distutils.core import setup
    from scipy_distutils.misc_util import merge_config_dicts
    from special_version import special_version

    config_dict = merge_config_dicts([configuration()] + \
                                     map(get_package_config,extra_packages))

    setup(version=special_version,
          **config_dict)

