#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.system_info import get_info, dict_append
    package = 'optimize'
    config = default_config_dict(package, parent_package)
    local_path = get_path(__name__,parent_path)    

    numpy_info = get_info('numpy',notfound_action=2)

    minpack = ('minpack',{'sources':
                          glob(os.path.join(local_path,'minpack','*.f'))})

    sources = ['_minpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext_args = {}
    dict_append(ext_args,
                name = dot_join(parent_package, package, '_minpack'),
                sources = sources, libraries = [minpack])
    dict_append(ext_args,**numpy_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    rootfind = glob(os.path.join(local_path,'Zeros','*.c'))
    roothead = os.path.join(local_path,'zeros.h')
    config['libraries'].append(('rootfind',{'sources':rootfind,
                                            'headers':roothead}))
    sources = ['zeros.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package, package, '_zeros'),
                sources=sources, libraries=['rootfind'])
    dict_append(ext_args,**numpy_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    lapack = get_info('lapack_opt')
    sources = ['lbfgsb.pyf','routines.f']
    sources = [os.path.join(local_path,'lbfgsb-0.9',x) for x in sources]
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package, package, "_lbfgsb"),
                sources=sources)
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**lapack)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['moduleTNC.c', 'tnc.c']
    sources = [os.path.join(local_path,'tnc',x) for x in sources]
    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package, package, "moduleTNC"),
                sources=sources)
    dict_append(ext_args,**numpy_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['cobyla.pyf','cobyla2.f','trstlp.f']
    sources = [os.path.join(local_path,'cobyla',x) for x in sources]
    ext = Extension(dot_join(parent_package, package, '_cobyla'), sources=sources)
    config['ext_modules'].append(ext)

    sources = ['minpack2.pyf', 'dcsrch.f', 'dcstep.f']
    sources = [os.path.join(local_path,'minpack2',x) for x in sources]

    ext_args = {}
    dict_append(ext_args,
                name=dot_join(parent_package, package, "minpack2"),
                sources=sources)
    dict_append(ext_args,**numpy_info)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
