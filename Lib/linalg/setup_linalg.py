#!/usr/bin/env python

from __future__ import nested_scopes
import os
import sys
import re
from distutils.dep_util import newer_group, newer
from glob import glob
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

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import dot_join, get_path, default_config_dict
    from scipy_distutils.system_info import get_info, dict_append, NotFoundError

    from interface_gen import generate_interface

    package = 'linalg'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    def local_join(*paths):
        return os.path.join(*((local_path,)+paths))

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError,'no lapack/blas resources found'

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0] 
    if atlas_version:
        print 'ATLAS version',atlas_version

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
    elif atlas_version>'3.4.0' and atlas_version<='3.5.12':
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
    ext_args = {'name': dot_join(parent_package,package,'fblas'),
                'sources': [generate_pyf,
                            local_join('src','fblaswrap.f')],
                'depends': map(local_join,['generic_fblas.pyf',
                                           'generic_fblas1.pyf',
                                           'generic_fblas2.pyf',
                                           'generic_fblas3.pyf',
                                           'interface_gen.py'])
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # cblas:
    ext_args = {'name': dot_join(parent_package,package,'cblas'),
                'sources': [generate_pyf],
                'depends': map(local_join,['generic_cblas.pyf',
                                           'generic_cblas1.pyf',
                                           'interface_gen.py'])
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # flapack:
    ext_args = {'name': dot_join(parent_package,package,'flapack'),
                'sources': [generate_pyf],
                'depends': map(local_join,['generic_flapack.pyf',
                                           'flapack_user_routines.pyf',
                                           'interface_gen.py'])
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # clapack:
    ext_args = {'name': dot_join(parent_package,package,'clapack'),
                'sources': [generate_pyf],
                'depends': map(local_join,['generic_clapack.pyf',
                                           'interface_gen.py'])
                }
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # _flinalg:
    ext_args = {'name':dot_join(parent_package,package,'_flinalg'),
                'sources':[local_join('src','det.f'),
                           local_join('src','lu.f')]
                }
    dict_append(ext_args,**lapack_opt)
    config['ext_modules'].append(Extension(**ext_args))

    # calc_lwork:
    ext_args = {'name':dot_join(parent_package,package,'calc_lwork'),
                'sources':[local_join('src','calc_lwork.f')],
                }
    dict_append(ext_args,**lapack_opt)
    config['ext_modules'].append(Extension(**ext_args))

    # atlas_version:
    ext_args = {'name':dot_join(parent_package,package,'atlas_version'),
                'sources':[os.path.join(local_path,'atlas_version.c')]}
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # iterative methods
    methods = ['BiCGREVCOM.f.src',
               'BiCGSTABREVCOM.f.src',
               'CGREVCOM.f.src',
               'CGSREVCOM.f.src',
#               'ChebyREVCOM.f.src',
               'GMRESREVCOM.f.src',
#               'JacobiREVCOM.f.src',
               'QMRREVCOM.f.src',
#               'SORREVCOM.f.src'
               ]
    Util = ['STOPTEST2.f.src','getbreak.f.src']
    sources = Util + methods + ['_iterative.pyf.src']
    ext_args = {
        'name': dot_join(parent_package,package,'_iterative'),
        'sources': [local_join('iterative',x) for x in sources]
        }
    dict_append(ext_args, **lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)    

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(parent_path=''))
