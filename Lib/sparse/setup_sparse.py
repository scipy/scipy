#!/usr/bin/env python
from __future__ import nested_scopes

import os
import sys
import glob

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
    from scipy_distutils.system_info import get_info,dict_append,NotFoundError

    package = 'sparse'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    def local_join(*paths):
        return os.path.join(*((local_path,)+paths))
    def local_glob(*names):
        return glob.glob(os.path.join(*((local_path,)+names)))

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError,'no lapack/blas resources found'

    if sys.platform=='win32':
        superlu_defs = [('NO_TIMER',1)]
    else:
        superlu_defs = []
    superlu_defs.append(('USE_VENDOR_BLAS',1))
    superlu = ('superlu_src',{'sources':local_glob('SuperLU','SRC','*.c'),'macros':superlu_defs})

    sparsekit = ('sparsekit_src',{'sources':local_glob('sparsekit','*.f')})

    #SuperLU/SRC/util.h  has been modifed to use these by default
    #macs = [('USER_ABORT','superlu_python_module_abort'),
    #        ('USER_MALLOC','superlu_python_module_malloc'),
    #        ('USER_FREE','superlu_python_module_free')]
    
    # Extension
    sources = ['_zsuperlumodule.c','_superlu_utils.c','_superluobject.c']
    ext_args = {'name':dot_join(parent_package,package,'_zsuperlu'),
                'sources':map(local_join,sources),
                'libraries': [superlu],
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['_dsuperlumodule.c','_superlu_utils.c','_superluobject.c']
    ext_args = {'name':dot_join(parent_package,package,'_dsuperlu'),
                'sources':map(local_join,sources),
                'libraries': [superlu],
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['_csuperlumodule.c','_superlu_utils.c','_superluobject.c']
    ext_args = {'name':dot_join(parent_package,package,'_csuperlu'),
                'sources':map(local_join,sources),
                'libraries': [superlu],
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    sources = ['_ssuperlumodule.c','_superlu_utils.c','_superluobject.c']
    ext_args = {'name':dot_join(parent_package,package,'_ssuperlu'),
                'sources':map(local_join,sources),
                'libraries': [superlu],
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

#    ext_args = {'name':dot_join(parent_package,package,'_sparsekit'),
#                'sources':[local_join('_sparsekit.pyf')],
#                'libraries' : [sparsekit]
#                }
#    dict_append(ext_args,**lapack_opt)
#    ext = Extension(**ext_args)
#    config['ext_modules'].append(ext)

    sources = ['spblas.f.src','spconv.f.src','sparsetools.pyf.src']
    sources = [local_join('sparsetools',x) for x in sources]
    ext = Extension(dot_join(parent_package, package, 'sparsetools'), sources=sources)
    config['ext_modules'].append(ext)
    
    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
