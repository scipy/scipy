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

tmpl_empty_clapack_pyf = '''
python module cblas
  usercode void empty_module(void) {}
  interface
    subroutine empty_module()
      intent(c) empty_module
    end subroutine empty_module
  end interface
end python module cblas
'''

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import dot_join, get_path, default_config_dict
    from scipy_distutils.system_info import get_info, dict_append

    package = 'blas'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    def local_join(*paths):
        return os.path.join(*((local_path,)+paths))
    def local_glob(path):
        return glob(os.path.join(local_path,path))

    numpy_info = get_info('numpy',notfound_action=2)
    blas_opt = get_info('blas_opt',notfound_action=2)

    atlas_version = ([v[3:-3] for k,v in blas_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0] 
    if atlas_version:
        print 'ATLAS version',atlas_version

    target_dir = ''
    skip_names = {'cblas':[],'fblas':[]}
    if skip_single_routines:
        target_dir = 'dbl'
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

    # fblas:
    ext_args = {}
    dict_append(ext_args,
                name = dot_join(parent_package,package,'fblas'),
                sources = [local_join('fblas.pyf.src'),
                           local_join('fblaswrap.f.src')],
                depends = [__file__]+local_glob('fblas_l?.pyf.src'),
                f2py_options = ['skip:']+skip_names['fblas']+[':'])
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**blas_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # cblas:
    def get_cblas_source(ext, build_dir):
        name = ext.name.split('.')[-1]
        assert name=='cblas',`name`
        if atlas_version is None:
            target = join(build_dir,target_dir,'cblas.pyf')
            from distutils.dep_util import newer
            if newer(__file__,target):
                f = open(source,'w')
                f.write(tmpl_empty_cblas_pyf)
                f.close()
        else:
            target = ext.depends[0]
            assert os.path.basename(target)=='cblas.pyf.src'
        return target

    ext_args = {}
    dict_append(ext_args,
                name = dot_join(parent_package,package,'cblas'),
                sources = [get_cblas_source],
                depends =  [local_join('cblas.pyf.src')] \
                + local_glob('cblas_l?.pyf.src'),
                f2py_options = ['skip:']+skip_names['cblas']+[':'])
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**blas_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration(parent_path=''))
