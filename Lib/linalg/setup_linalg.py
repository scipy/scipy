#!/usr/bin/env python

import os,sys
from distutils import dep_util
from glob import glob
import warnings

# XXX: system_info.py should set the following values.

# Changing atlas_version_pre_3_3 value requires
#  rm -f clapack.pyf
# before rebuilding.
atlas_version_pre_3_3 = 0

# Changing skip_single_routines value requires
#  rm -f {clapack,flapack,cblas,fblas}.pyf
# before rebuilding.
skip_single_routines = 0

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
     LapackSrcNotFoundError,BlasSrcNotFoundError

def configuration(parent_package=''):
    package = 'linalg'
    from interface_gen import generate_interface
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__)

    atlas_info = get_info('atlas')
    #atlas_info = {} # uncomment if ATLAS is available but want to use
                     # Fortran LAPACK/BLAS; useful for testing
    f_libs = []
    blas_info,lapack_info = {},{}
    if not atlas_info:
        warnings.warn(AtlasNotFoundError.__doc__)
        blas_info = get_info('blas')
        #blas_info = {} # test building BLAS from sources.
        if not blas_info:
            warnings.warn(BlasNotFoundError.__doc__)
            blas_src_info = get_info('blas_src')
            if not blas_src_info:
                raise BlasSrcNotFoundError,BlasSrcNotFoundError.__doc__
            dict_append(blas_info,libraries=['blas_src'])
            f_libs.append(fortran_library_item(\
                'blas_src',blas_src_info['sources'],
                ))
        lapack_info = get_info('lapack')
        #lapack_info = {} # test building LAPACK from sources.
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

    skip_names = {'clapack':[],'flapack':[],'cblas':[],'fblas':[]}
    if skip_single_routines:
        skip_names['clapack'].extend(\
            'sgesv cgesv sgetrf cgetrf sgetrs cgetrs sgetri cgetri'\
            ' sposv cposv spotrf cpotrf spotrs cpotrs spotri cpotri'\
            ' slauum clauum strtri ctrtri'.split())
        skip_names['flapack'].extend(skip_names['clapack'])
        skip_names['flapack'].extend(\
            'sgesdd cgesdd sgelss cgelss sgeqrf cgeqrf sgeev cgeev'\
            ' sgegv cgegv ssyev cheev slaswp claswp sgees cgees'
            ' sggev cggev'.split())
        skip_names['cblas'].extend('saxpy caxpy'.split())
        skip_names['fblas'].extend(skip_names['cblas'])
        skip_names['fblas'].extend(\
            'srotg crotg srotmg srot csrot srotm sswap cswap sscal cscal'\
            ' csscal scopy ccopy sdot cdotu cdotc snrm2 scnrm2 sasum scasum'\
            ' isamax icamax sgemv cgemv chemv ssymv strmv ctrmv'\
            ' sgemm cgemm'.split())
    if atlas_version_pre_3_3:
        skip_names['clapack'].extend(\
            'sgetri dgetri cgetri zgetri spotri dpotri cpotri zpotri'\
            ' slauum dlauum clauum zlauum strtri dtrtri ctrtri ztrtri'.split())

    for mod_name,sources in mod_sources.items():
        sources = [os.path.join(local_path,s) for s in sources]
        mod_file = os.path.join(local_path,mod_name+'.pyf')
        if dep_util.newer_group(sources,mod_file):
            generate_interface(mod_name,sources[0],mod_file,
                               skip_names.get(mod_name,[]))
        sources = filter(lambda s:s[-4:]!='.pyf',sources)
        ext_args = {'name':dot_join(parent_package,package,mod_name),
                    'sources':[mod_file]+sources}
        dict_append(ext_args,**atlas_info)
        ext = Extension(**ext_args)
        ext.need_fcompiler_opts = 1
        config['ext_modules'].append(ext)

    flinalg = []
    for f in ['det.f','lu.f', #'wrappers.c','inv.f',
              ]:
        flinalg.append(os.path.join(local_path,'src',f))
    ext_args = {'name':dot_join(parent_package,package,'_flinalg'),
                'sources':flinalg}
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    ext_args = {'name':dot_join(parent_package,package,'calc_lwork'),
                'sources':[os.path.join(local_path,'src','calc_lwork.f')],
                }
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    config['fortran_libraries'].extend(f_libs)
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
