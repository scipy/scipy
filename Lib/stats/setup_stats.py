#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package='',parent_path=None):
    package = 'stats'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package, parent_package)

    statlib = glob(os.path.join(local_path, 'statlib','*.f'))

    config['fortran_libraries'].append(('statlib',{'sources':statlib}))
    
    # Extension
    sources = ['randmodule.c','ranlib_all.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'rand'),sources)
    config['ext_modules'].append(ext)

    ext_args = {'name':dot_join(parent_package,package,'statlib'),
                'sources':[os.path.join(local_path,'statlib.pyf')],
                'f2py_options':['--no-wrap-functions'],
                #'define_macros':[('F2PY_REPORT_ATEXIT_DISABLE',None)],
                'libraries' : ['statlib']
                }
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    # add futil module
    sources = ['futil.f']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'futil'),sources)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
