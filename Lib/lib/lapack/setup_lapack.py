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


def configuration(parent_package='',parent_path=None):
    from scipy_distutils.core import Extension
    from scipy_distutils.misc_util import dot_join, get_path, default_config_dict
    from scipy_distutils.system_info import get_info, dict_append

    package = 'lapack'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)
    def local_join(*paths):
        return os.path.join(*((local_path,)+paths))
    def local_glob(path):
        return glob(os.path.join(local_path,path))

    numpy_info = get_info('numpy',notfound_action=2)
    lapack_opt = get_info('lapack_opt',notfound_action=2)

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print 'ATLAS version',atlas_version

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
    elif atlas_version>'3.4.0' and atlas_version<='3.5.12':
        skip_names['clapack'].extend('cpotrf zpotrf'.split())

    # flapack:
    ext_args = {}
    dict_append(ext_args,
                name = dot_join(parent_package,package,'flapack'),
                sources = [local_join('flapack.pyf.src')],
                depends = [__file__]+local_glob('flapack_*.pyf.src'),
                f2py_options = ['skip:']+skip_names['flapack']+[':'])
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)


    # clapack:
    def get_clapack_source(ext, build_dir):
        name = ext.name.split('.')[-1]
        assert name=='clapack',`name`
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

    ext_args = {}
    dict_append(ext_args,
                name = dot_join(parent_package,package,'clapack'),
                sources = [get_clapack_source],
                depends =  [local_join('clapack.pyf.src')] \
                + local_glob('clapack_*.pyf.src'),
                f2py_options = ['skip:']+skip_names['clapack']+[':'])
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    # calc_lwork:
    ext_args = {'name':dot_join(parent_package,package,'calc_lwork'),
                'sources':[local_join('calc_lwork.f')],
                }
    dict_append(ext_args,**numpy_info)
    dict_append(ext_args,**lapack_opt)
    config['ext_modules'].append(Extension(**ext_args))

    # atlas_version:
    ext_args = {'name':dot_join(parent_package,package,'atlas_version'),
                'sources':[os.path.join(local_path,'atlas_version.c')]}
    dict_append(ext_args,**lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup

    setup(**configuration(parent_path=''))
