#!/usr/bin/env python

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils import system_info


def configuration(parent_package='',parent_path=None):
    package = 'optimize'
    config = default_config_dict(package, parent_package)
    local_path = get_path(__name__,parent_path)    

    minpack = ('minpack',{'sources':
                          glob(os.path.join(local_path,'minpack','*.f'))})

    sources = ['_minpackmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package, package, '_minpack'),
                    sources, libraries = [minpack])
    config['ext_modules'].append(ext)

    rootfind = glob(os.path.join(local_path,'Zeros','*.c'))
    roothead = os.path.join(local_path,'zeros.h')
    config['libraries'].append(('rootfind',{'sources':rootfind,
                                            'headers':roothead}))
    sources = ['zeros.c']
    sources = [os.path.join(local_path,x) for x in sources]    
    ext = Extension(dot_join(parent_package, package, '_zeros'),
                    sources, libraries=['rootfind'])
    config['ext_modules'].append(ext)

    lapack = system_info.get_info('lapack_opt')
    sources = ['lbfgsb.pyf','routines.f']
    sources = [os.path.join(local_path,'lbfgsb-0.9',x) for x in sources]
    ext = Extension(dot_join(parent_package, package, "_lbfgsb"),
		    sources=sources, **lapack)
    config['ext_modules'].append(ext)

    sources = ['moduleTNC.c', 'tnc.c']
    sources = [os.path.join(local_path,'tnc',x) for x in sources]
    ext = Extension(dot_join(parent_package, package, 'moduleTNC'), 
		    sources=sources)
    config['ext_modules'].append(ext)

    sources = ['cobyla.pyf','cobyla2.f','trstlp.f']
    sources = [os.path.join(local_path,'cobyla',x) for x in sources]
    ext = Extension(dot_join(parent_package, package, '_cobyla'), sources=sources)
    config['ext_modules'].append(ext)

    sources = ['minpack2.pyf', 'dcsrch.f', 'dcstep.f']
    sources = [os.path.join(local_path,'minpack2',x) for x in sources]    
    ext = Extension(dot_join(parent_package, package, 'minpack2'), sources=sources)
    config['ext_modules'].append(ext)
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
