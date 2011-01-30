#!/usr/bin/env python

import os
from glob import glob

#-------------------
# To skip wrapping single precision atlas/lapack routines, set
# the following flag to True:

skip_single_routines = 0

#--------------------

tmpl_empty_clapack_pyf = '''
python module clapack
  usercode void empty_module(void) {}
  interface
    subroutine empty_module()
      intent(c) empty_module
    end subroutine empty_module
  end interface
end python module clapack
'''


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info

    config = Configuration('lapack',parent_package,top_path)

    lapack_opt = get_info('lapack_opt',notfound_action=2)

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print ('ATLAS version: %s' % atlas_version)

    target_dir = ''
    skip_names = {'clapack':[],'flapack':[]}
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

    if atlas_version=='3.2.1_pre3.3.6':
        target_dir = os.path.join(target_dir,'atlas321')
        skip_names['clapack'].extend(\
            'sgetri dgetri cgetri zgetri spotri dpotri cpotri zpotri'\
            ' slauum dlauum clauum zlauum strtri dtrtri ctrtri ztrtri'.split())
    elif atlas_version and atlas_version>'3.4.0' and atlas_version<='3.5.12':
        skip_names['clapack'].extend('cpotrf zpotrf'.split())

    # flapack:
    config.add_extension('flapack',
                         sources = ['flapack.pyf.src'],
                         depends = [__file__,'flapack_*.pyf.src'],
                         f2py_options = ['skip:']+skip_names['flapack']+[':'],
                         extra_info = lapack_opt
                         )

    # clapack:
    def get_clapack_source(ext, build_dir):
        name = ext.name.split('.')[-1]
        assert name=='clapack', repr(name)
        if atlas_version is None:
            target = os.path.join(build_dir,target_dir,'clapack.pyf')
            from distutils.dep_util import newer
            if newer(__file__,target):
                f = open(target,'w')
                f.write(tmpl_empty_clapack_pyf)
                f.close()
        else:
            target = ext.depends[0]
            assert os.path.basename(target)=='clapack.pyf.src'
        return target

    config.add_extension('clapack',
                         sources = [get_clapack_source],
                         depends = ['clapack.pyf.src'],
                         f2py_options = ['skip:']+skip_names['clapack']+[':'],
                         extra_info = lapack_opt
                         )

    # calc_lwork:
    config.add_extension('calc_lwork',
                         sources = ['calc_lwork.f'],
                         extra_info = lapack_opt
                         )

    # atlas_version:
    if os.name == 'nt' and 'FPATH' in os.environ:
        define_macros = [('NO_ATLAS_INFO', 1)]
    else:
        define_macros = []

    config.add_extension('atlas_version',
                         sources = ['atlas_version.c'],
                         extra_info = lapack_opt,
                         define_macros = define_macros
                         )

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
