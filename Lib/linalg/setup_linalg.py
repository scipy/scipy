#!/usr/bin/env python

import os,sys
from distutils import dep_util
from glob import glob
import warnings

if sys.platform == 'win32':
    # force g77 for now
    # XXX: g77 is forced already in scipy/setup.py
    #      So, is this redundant code in what follows?
    from scipy_distutils.mingw32_support import *
    from scipy_distutils.command import build_flib
    build_flib.all_compilers = [build_flib.gnu_fortran_compiler]

from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils.misc_util import fortran_library_item

from scipy_distutils.system_info import get_info,dict_append,\
     AtlasNotFoundError,LapackNotFoundError,BlasNotFoundError,\
     LapackSrcNotFoundError

def configuration(parent_package=''):
    from interface_gen import generate_interface
    config = default_config_dict('linalg',parent_package)
    local_path = get_path(__name__)
    test_path = os.path.join(local_path,'tests')

    config['packages'].append(dot_join(parent_package,'linalg.tests'))
    config['package_dir']['linalg.tests'] = test_path

    atlas_info = get_info('atlas')
    #atlas_info = {} # uncomment if ATLAS is available but want to use
                     # Fortran LAPACK/ATLAS; useful for testing
    f_libs = []
    blas_info,lapack_info,lapack_src_info = {},{},{}
    if not atlas_info:
        warnings.warn(AtlasNotFoundError.__doc__)
        blas_info = get_info('blas')
        lapack_info = get_info('lapack')
        #lapack_info = {} # test building LAPACK from sources.
        if not blas_info:
            #TODO: build from blas sources (see lapack_src below)
            raise BlasNotFoundError,BlasNotFoundError.__doc__
        if not lapack_info:
            warnings.warn(LapackNotFoundError.__doc__)
            lapack_src_info = get_info('lapack_src')
            if not lapack_src_info:
                raise LapackSrcNotFoundError,LapackSrcNotFoundError.__doc__
            dict_append(lapack_info,libraries=['lapack_src'])
            f_libs.append(fortran_library_item(\
                'lapack_src',lapack_src_info['sources'],
                ))

    mod_sources = {}
    if atlas_info or blas_info:
        mod_sources['fblas'] = ['generic_fblas.pyf',
                                'generic_fblas1.pyf',
                                'generic_fblas2.pyf',
                                'generic_fblas3.pyf',
                                os.path.join('src','fblaswrap.f'),
                                ]
    if atlas_info or lapack_info:
        mod_sources['flapack'] = ['generic_flapack.pyf']
    if atlas_info:
        mod_sources['cblas'] = ['generic_cblas.pyf',
                                'generic_cblas1.pyf']
        mod_sources['clapack'] = ['generic_clapack.pyf']
    else:
        dict_append(atlas_info,**lapack_info)
        dict_append(atlas_info,**blas_info)

    for mod_name,sources in mod_sources.items():
        sources = [os.path.join(local_path,s) for s in sources]
        mod_file = os.path.join(local_path,mod_name+'.pyf')
        if dep_util.newer_group(sources,mod_file):
            generate_interface(mod_name,sources[0],mod_file)
        sources = filter(lambda s:s[-4:]!='.pyf',sources)
        ext_args = {'name':dot_join(parent_package,'linalg',mod_name),
                    'sources':[mod_file]+sources}
        dict_append(ext_args,**atlas_info)
        ext = Extension(**ext_args)
        ext.need_fcompiler_opts = 1
        config['ext_modules'].append(ext)

    flinalg = []
    for f in ['det.f','lu.f', #'wrappers.c','inv.f',
              ]:
        flinalg.append(os.path.join(local_path,'src',f))
    ext_args = {'name':dot_join(parent_package,'linalg','_flinalg'),
                'sources':flinalg}
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    ext_args = {'name':dot_join(parent_package,'linalg','calc_lwork'),
                'sources':[os.path.join(local_path,'src','calc_lwork.f')],
                }
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    config['fortran_libraries'].extend(f_libs)
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
