#!/usr/bin/env python

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

def needs_cblas_wrapper(info):
    """Returns true if needs c wrapper around cblas for calling from
    fortran."""
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

tmpl_empty_cblas_pyf = '''
python module cblas
  usercode void empty_module(void) {}
  interface
    subroutine empty_module()
      intent(c) empty_module
    end subroutine empty_module
  end interface
end python module cblas
'''

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('blas',parent_package,top_path)

    blas_opt = get_info('blas_opt',notfound_action=2)

    atlas_version = ([v[3:-3] for k,v in blas_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print ('ATLAS version: %s' % atlas_version)

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

    depends = [__file__, 'fblas_l?.pyf.src', 'fblas.pyf.src','fblaswrap.f.src',
               'fblaswrap_veclib_c.c.src']
    # fblas:
    if needs_cblas_wrapper(blas_opt):
        sources = ['fblas.pyf.src', 'fblaswrap_veclib_c.c.src'],
    else:
        sources = ['fblas.pyf.src','fblaswrap.f.src']
    config.add_extension('fblas',
                         sources = sources,
                         depends = depends,
                         f2py_options = ['skip:']+skip_names['fblas']+[':'],
                         extra_info = blas_opt
                         )
    # cblas:
    def get_cblas_source(ext, build_dir):
        name = ext.name.split('.')[-1]
        assert name=='cblas', repr(name)
        if atlas_version is None:
            target = join(build_dir,target_dir,'cblas.pyf')
            from distutils.dep_util import newer
            if newer(__file__,target):
                f = open(target,'w')
                f.write(tmpl_empty_cblas_pyf)
                f.close()
        else:
            target = ext.depends[0]
            assert os.path.basename(target)=='cblas.pyf.src'
        return target

    config.add_extension('cblas',
                         sources = [get_cblas_source],
                         depends = ['cblas.pyf.src','cblas_l?.pyf.src'],
                         f2py_options = ['skip:']+skip_names['cblas']+[':'],
                         extra_info = blas_opt
                         )

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
