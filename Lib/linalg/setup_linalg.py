#!/usr/bin/env python

import os
import sys
import re
from distutils import dep_util,dir_util
from glob import glob
import warnings
import shutil

#-------------------
# To skip wrapping single precision atlas/lapack/blas routines, set
# the following flag to True:
skip_single_routines = 0

# Some OS distributions (e.g. Redhat, Suse) provide a blas library that
# is built using incomplete blas sources that come with lapack tar-ball.
# In order to use such a library in scipy.linalg, the following flag
# must be set to True:
using_lapack_blas = 0

#--------------------

if os.name == 'nt':
    def run_command(command):
        """ not sure how to get exit status on nt. """
        in_pipe,out_pipe = os.popen4(command)
        in_pipe.close()
        text = out_pipe.read()
        return 0, text
else:
    import commands
    run_command = commands.getstatusoutput

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import fortran_library_item, dot_join,\
         SourceGenerator, get_path, default_config_dict, get_build_temp
    from scipy_distutils.system_info import get_info,dict_append,\
         AtlasNotFoundError,LapackNotFoundError,BlasNotFoundError,\
         LapackSrcNotFoundError,BlasSrcNotFoundError

    package = 'linalg'
    from interface_gen import generate_interface
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    abs_local_path = os.path.abspath(local_path)

    atlas_info = get_info('atlas')
    #atlas_info = {} # uncomment if ATLAS is available but want to use
                     # Fortran LAPACK/BLAS; useful for testing

    f_libs = []
    atlas_version = None
    temp_path = os.path.join(get_build_temp(),'linalg','atlas_version')
    dir_util.mkpath(temp_path,verbose=1)
    atlas_version_file = os.path.join(temp_path,'atlas_version')

    if atlas_info:
        if os.path.isfile(atlas_version_file):
            atlas_version = open(atlas_version_file).read()
            print 'ATLAS version',atlas_version

    if atlas_info and atlas_version is None:
        # Try to determine ATLAS version
        shutil.copy(os.path.join(local_path,'atlas_version.c'),temp_path)
        cur_dir = os.getcwd()
        os.chdir(temp_path)
        cmd = '%s %s build_ext --inplace --force'%\
              (sys.executable,
               os.path.join(abs_local_path,'setup_atlas_version.py'))
        print cmd
        s,o=run_command(cmd)
        if not s:
            cmd = sys.executable+' -c "import atlas_version"'
            print cmd
            s,o=run_command(cmd)
            if not s:
                m = re.match(r'ATLAS version (?P<version>\d+[.]\d+[.]\d+)',o)
                if m:
                    atlas_version = m.group('version')
                    print 'ATLAS version',atlas_version
            if atlas_version is None:
                if re.search(r'undefined symbol: ATL_buildinfo',o,re.M):
                    atlas_version = '3.2.1_pre3.3.6'
                    print 'ATLAS version',atlas_version
                else:
                    print o
        else:
            print o
        os.chdir(cur_dir)
        if atlas_version is None:
            print 'Failed to determine ATLAS version'
        else:
            f = open(atlas_version_file,'w')
            f.write(atlas_version)
            f.close()

    if atlas_info:
        if ('ATLAS_WITHOUT_LAPACK',None) in atlas_info.get('define_macros',[]):
            lapack_info = get_info('lapack')
            if not lapack_info:
                warnings.warn(LapackNotFoundError.__doc__)
                lapack_src_info = get_info('lapack_src')
                if not lapack_src_info:
                    raise LapackSrcNotFoundError,LapackSrcNotFoundError.__doc__
                dict_append(lapack_info,libraries=['lapack_src'])
                f_libs.append(fortran_library_item(\
                'lapack_src',lapack_src_info['sources'],
                ))
            dict_append(lapack_info,**atlas_info)
            atlas_info = lapack_info

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
                'blas_src',blas_src_info['sources'] + \
                [os.path.join(local_path,'src','fblaswrap.f')],
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

        dict_append(atlas_info,**lapack_info)
        dict_append(atlas_info,**blas_info)

    target_dir = ''
    skip_names = {'clapack':[],'flapack':[],'cblas':[],'fblas':[]}
    if skip_single_routines:
        target_dir = 'dbl'
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
    if using_lapack_blas:
        target_dir = os.path.join(target_dir,'blas')
        skip_names['fblas'].extend(\
            'drotmg srotmg drotm srotm'.split())

    if atlas_version=='3.2.1_pre3.3.6':
        target_dir = os.path.join(target_dir,'atlas321')
        skip_names['clapack'].extend(\
            'sgetri dgetri cgetri zgetri spotri dpotri cpotri zpotri'\
            ' slauum dlauum clauum zlauum strtri dtrtri ctrtri ztrtri'.split())

    # atlas_version:
    ext_args = {'name':dot_join(parent_package,package,'atlas_version'),
                'sources':[os.path.join(local_path,'atlas_version.c')]}
    if atlas_info:
        ext_args['libraries'] = [atlas_info['libraries'][-1]]
        ext_args['library_dirs'] = atlas_info['library_dirs'][:]
        ext_args['define_macros'] = [('ATLAS_INFO','"%s"' % atlas_version)]
    else:
        ext_args['define_macros'] = [('NO_ATLAS_INFO',1)]
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # In case any of atlas|lapack|blas libraries are not available
    def generate_empty_pyf(target,sources,generator,skips):
        name = os.path.basename(target)[:-4]
        f = open(target,'w')
        f.write('python module '+name+'\n')
        f.write('usercode void empty_module(void) {}\n')
        f.write('interface\n')
        f.write('subroutine empty_module()\n')
        f.write('intent(c) empty_module\n')
        f.write('end subroutine empty_module\n')
        f.write('end interface\nend python module'+name+'\n')
        f.close()

    # fblas:
    def generate_fblas_pyf(target,sources,generator,skips):
        generator('fblas',sources[0],target,skips)
    if not (blas_info or atlas_info):
        generate_fblas_pyf = generate_empty_pyf
    sources = ['generic_fblas.pyf',
               'generic_fblas1.pyf',
               'generic_fblas2.pyf',
               'generic_fblas3.pyf',
               os.path.join('src','fblaswrap.f')]
    sources = [os.path.join(local_path,s) for s in sources]
    fblas_pyf = SourceGenerator(generate_fblas_pyf,
                                os.path.join(target_dir,'fblas.pyf'),
                                sources,generate_interface,
                                skip_names['fblas'])
    ext_args = {'name':dot_join(parent_package,package,'fblas'),
                'sources':[fblas_pyf,sources[-1]]}
    dict_append(ext_args,**atlas_info)
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    # cblas:
    def generate_cblas_pyf(target,sources,generator,skips):
        generator('cblas',sources[0],target,skips)
    if not atlas_info:
        generate_cblas_pyf = generate_empty_pyf
    sources = ['generic_cblas.pyf',
               'generic_cblas1.pyf']
    sources = [os.path.join(local_path,s) for s in sources]
    cblas_pyf = SourceGenerator(generate_cblas_pyf,
                                os.path.join(target_dir,'cblas.pyf'),
                                sources,generate_interface,
                                skip_names['cblas'])
    ext_args = {'name':dot_join(parent_package,package,'cblas'),
                'sources':[cblas_pyf]}
    dict_append(ext_args,**atlas_info)
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    # flapack:
    def generate_flapack_pyf(target,sources,generator,skips):
        generator('flapack',sources[0],target,skips)
    if not (lapack_info or atlas_info):
        generate_flapack_pyf = generate_empty_pyf
    sources = ['generic_flapack.pyf','flapack_user_routines.pyf']
    sources = [os.path.join(local_path,s) for s in sources]
    flapack_pyf = SourceGenerator(generate_flapack_pyf,
                                  os.path.join(target_dir,'flapack.pyf'),
                                  sources,generate_interface,
                                  skip_names['flapack'])
    ext_args = {'name':dot_join(parent_package,package,'flapack'),
                'sources':[flapack_pyf]}
    dict_append(ext_args,**atlas_info)
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    # clapack:
    def generate_clapack_pyf(target,sources,generator,skips):
        generator('clapack',sources[0],target,skips)
    if not atlas_info:
        generate_cblas_pyf = generate_empty_pyf
    sources = ['generic_clapack.pyf']
    sources = [os.path.join(local_path,s) for s in sources]
    clapack_pyf = SourceGenerator(generate_clapack_pyf,
                                  os.path.join(target_dir,'clapack.pyf'),
                                  sources,generate_interface,
                                  skip_names['clapack'])
    ext_args = {'name':dot_join(parent_package,package,'clapack'),
                'sources':[clapack_pyf]}
    dict_append(ext_args,**atlas_info)
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    config['ext_modules'].append(ext)

    # _flinalg:
    flinalg = []
    for f in ['det.f','lu.f', #'wrappers.c','inv.f',
              ]:
        flinalg.append(os.path.join(local_path,'src',f))
    ext_args = {'name':dot_join(parent_package,package,'_flinalg'),
                'sources':flinalg}
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    # calc_lwork:
    ext_args = {'name':dot_join(parent_package,package,'calc_lwork'),
                'sources':[os.path.join(local_path,'src','calc_lwork.f')],
                }
    dict_append(ext_args,**atlas_info)
    config['ext_modules'].append(Extension(**ext_args))

    config['fortran_libraries'].extend(f_libs)
    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration())
