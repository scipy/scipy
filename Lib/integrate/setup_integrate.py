#!/usr/bin/env python

from __future__ import nested_scopes
import os
import glob

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
    from scipy_distutils.system_info import get_info,dict_append,NotFoundError

    package = 'integrate'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    def local_join(*paths):
        return os.path.join(*((local_path,)+paths))
    def local_glob(*names):
        return glob.glob(os.path.join(*((local_path,)+names)))

    blas_opt = get_info('blas_opt')
    if not blas_opt:
        raise NotFoundError,'no blas resources found'

    linpack_lite = ('linpack_lite_src',{'sources':local_glob('linpack_lite','*.f')})
    mach = ('mach_src',{'sources':local_glob('mach','*.f')})
    quadpack = ('quadpack_src',{'sources':local_glob('quadpack','*.f')})
    odepack = ('odepack_src',{'sources':local_glob('odepack','*.f')})

    # should we try to weed through files and replace with calls to
    # LAPACK routines?
    # Yes, someday...

    # Extensions
    # quadpack:
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package,package,'_quadpack'),
                sources = [local_join('_quadpackmodule.c')],
                libraries = [quadpack,linpack_lite,mach],
                )
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # odepack
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package,package,'_odepack'),
                sources = [local_join('_odepackmodule.c')],
                libraries = [odepack,linpack_lite,mach],
                )
    dict_append(ext_args,**blas_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # vode
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package,package,'vode'),
                sources = [local_join('vode.pyf')],
                libraries = [odepack,linpack_lite,mach],
                )
    dict_append(ext_args,**blas_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
