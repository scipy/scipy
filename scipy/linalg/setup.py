#!/usr/bin/env python

import os
from distutils.dep_util import newer_group, newer
from os.path import join

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

def needs_cblas_wrapper(info):
    """Returns true if needs c wrapper around cblas for calling from
    fortran."""
    import re
    r_accel = re.compile("Accelerate")
    r_vec = re.compile("vecLib")
    res = False
    try:
        tmpstr = info['extra_link_args']
        for i in tmpstr:
            if r_accel.search(i) or r_vec.search(i):
                res = True
    except KeyError:
        pass

    return res

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError

    from numpy.distutils.misc_util import Configuration

    from interface_gen import generate_interface

    config = Configuration('linalg',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print ('ATLAS version: %s' % atlas_version)

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
        target_dir = join(target_dir,'blas')
        skip_names['fblas'].extend(\
            'drotmg srotmg drotm srotm'.split())

    if atlas_version=='3.2.1_pre3.3.6':
        target_dir = join(target_dir,'atlas321')
        skip_names['clapack'].extend(\
            'sgetri dgetri cgetri zgetri spotri dpotri cpotri zpotri'\
            ' slauum dlauum clauum zlauum strtri dtrtri ctrtri ztrtri'.split())
    elif atlas_version and atlas_version>'3.4.0' and atlas_version<='3.5.12':
        skip_names['clapack'].extend('cpotrf zpotrf'.split())

    def generate_pyf(extension, build_dir):
        name = extension.name.split('.')[-1]
        target = join(build_dir,target_dir,name+'.pyf')
        if name[0]=='c' and atlas_version is None and newer(__file__,target):
            f = open(target,'w')
            f.write('python module '+name+'\n')
            f.write('usercode void empty_module(void) {}\n')
            f.write('interface\n')
            f.write('subroutine empty_module()\n')
            f.write('intent(c) empty_module\n')
            f.write('end subroutine empty_module\n')
            f.write('end interface\nend python module'+name+'\n')
            f.close()
            return target
        if newer_group(extension.depends,target):
            generate_interface(name,
                               extension.depends[0],
                               target,
                               skip_names[name])
        return target


    # fblas:
    if needs_cblas_wrapper(lapack_opt):
        sources = ['fblas.pyf.src', join('src', 'fblaswrap_veclib_c.c')],
    else:
        sources = ['fblas.pyf.src', join('src', 'fblaswrap.f')]

    # Note: `depends` needs to include fblaswrap(_veclib) for both files to be
    # included by "python setup.py sdist"
    config.add_extension('fblas',
                         sources = sources,
                         depends = ['fblas_l?.pyf.src',
                                    join('src', 'fblaswrap_veclib_c.c'),
                                    join('src', 'fblaswrap.f')],
                         extra_info = lapack_opt
                         )

    # cblas:
    config.add_extension('cblas',
                         sources = [generate_pyf],
                         depends = ['generic_cblas.pyf',
                                    'generic_cblas1.pyf',
                                    'interface_gen.py'],
                         extra_info = lapack_opt
                         )

    # flapack:
    config.add_extension('flapack',
                         sources = [generate_pyf],
                         depends = ['generic_flapack.pyf',
                                    'flapack_user_routines.pyf',
                                    'interface_gen.py'],
                         extra_info = lapack_opt
                         )

    # clapack:
    config.add_extension('clapack',
                         sources = [generate_pyf],
                         depends = ['generic_clapack.pyf',
                                    'interface_gen.py'],
                         extra_info = lapack_opt
                         )

    # _flinalg:
    config.add_extension('_flinalg',
                         sources = [join('src','det.f'),join('src','lu.f')],
                         extra_info = lapack_opt
                         )

    # calc_lwork:
    config.add_extension('calc_lwork',
                         [join('src','calc_lwork.f')],
                         extra_info = lapack_opt
                         )

    # atlas_version:
    if os.name == 'nt' and 'FPATH' in os.environ:
        define_macros = [('NO_ATLAS_INFO', 1)]
    else:
        define_macros = []

    config.add_extension('atlas_version',
                         ['atlas_version.c'],
                         extra_info = lapack_opt,
                         define_macros = define_macros
                         )

    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(top_path='').todict())
