#!/usr/bin/env python

import os, sys
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils.system_info import dict_append
import shutil

def configuration(parent_package=''):
    package = 'special'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__)

    c_misc = glob(os.path.join(local_path,'c_misc','*.c'))
    cephes = glob(os.path.join(local_path,'cephes','*.c'))
    mach = glob(os.path.join(local_path,'mach','*.f'))
    amos = glob(os.path.join(local_path,'amos','*.f'))
    toms = glob(os.path.join(local_path,'toms','*.f'))
    cdf = glob(os.path.join(local_path,'cdflib','*.f'))
    specfun = glob(os.path.join(local_path, 'specfun','*.f'))
    
    # C libraries
    config['libraries'].append(('c_misc',{'sources':c_misc}))
    config['libraries'].append(('cephes',{'sources':cephes}))
    
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
                                 'cdf', 'specfun']
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

    # Test to see if big or little-endian machine and get correct default
    #   mconf.h module.
    cephes_path = os.path.join(local_path, 'cephes')
    if sys.byteorder == "little":
        print "### Little Endian detected ####"
        shutil.copy2(os.path.join(cephes_path,'mconf_LE.h'),os.path.join(cephes_path,'mconf.h'))
    else:
        print "### Big Endian detected ####"
        shutil.copy2(os.path.join(cephes_path,'mconf_BE.h'),os.path.join(cephes_path,'mconf.h'))

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
