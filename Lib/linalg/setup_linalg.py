#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os, sys, string
# this isn't in scipy_distutils -- need to figure out clean way to do that
from distutils import dep_util
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
import interface_gen
reload(interface_gen)
from scipy_distutils.atlas_info import get_atlas_info
from scipy_distutils.system_info import dict_append,AtlasNotFoundError

# needed now for pyf_extensions

if os.name == 'nt':
    from scipy_distutils.mingw32_support import *

def configuration(parent_package=''):
    config = default_config_dict('linalg',parent_package)
    local_path = get_path(__name__)
    test_path = os.path.join(local_path,'tests')

    config['packages'].append(dot_join(parent_package,'linalg.tests'))
    config['package_dir']['linalg.tests'] = test_path
    config['ext_modules'].extend(extensions(dot_join(parent_package,'linalg')))
    return config
       
def generic_extension(mod_name,sources,parent_package='',use_underscore = 1):
    atlas_info = get_atlas_info()
    if not atlas_info:
        raise AtlasNotFoundError,AtlasNotFoundError.__doc__
    
    local_path = get_path(__name__)
    mod_file = mod_name + '.pyf'

    #if parent_package and parent_package[-1] != '.':
    #    parent_package = parent_package + '.'    
    sources = [os.path.join(local_path,x) for x in sources]

    from distutils.dir_util import mkpath
    output_path = os.path.join('build','generated_pyfs')
    mkpath(output_path)
    mod_file = os.path.join(output_path,mod_file)

    # choose the
    gen_function = eval('interface_gen.generate_'+mod_name)
    
    if use_underscore:
        define_macros = [('F2PY_REPORT_ATEXIT_DISBLE',None)]
        define_macros = [('F2PY_REPORT_ATEXIT_DISABLE',None)]
    else:
        define_macros=[('NO_APPEND_FORTRAN',1),
                       ('F2PY_REPORT_ATEXIT_DISBLE',None)]
	define_macros=[('NO_APPEND_FORTRAN',1),
	               ('F2PY_REPORT_ATEXIT_DISABLE',None)]
            
    if dep_util.newer_group(sources,mod_file):
        gen_function(local_path,output_path)

    print 'file:', mod_file
    ext_args = {'name':dot_join(parent_package,mod_name),
                'sources':[mod_file],
                'f2py_options':['--no-wrap-functions'],
                'define_macros':define_macros,
                }
    dict_append(ext_args,**atlas_info)
    ext = Extension(**ext_args)
    ext.need_fcompiler_opts = 1
    return ext
    ext = Extension(dot_join(parent_package,mod_name),[mod_file],
                    library_dirs=atlas_library_dirs,
                    libraries = lapack_libraries,
                    define_macros=define_macros,
                    #XXX --no-wrap-functions is only needed for clapack
                    f2py_options=['--no-wrap-functions'],
                    )

    return ext                

def extensions(parent_package = ''):
    ext_modules = []
    
    # fortran lapack interface
    mod = 'flapack'
    sources = ['generic_lapack.pyf']
    ext = generic_extension(mod,sources,parent_package,use_underscore=1)
    ext_modules.append(ext)
    
    # C lapack interface
    mod = 'clapack'
    sources = ['generic_lapack.pyf']
    ext = generic_extension(mod,sources,parent_package,use_underscore=0)
    ext_modules.append(ext)

    # fblas lapack interface
    mod = 'fblas'
    sources = ['generic_blas1.pyf','generic_blas2.pyf','generic_blas3.pyf']
    ext = generic_extension(mod,sources,parent_package,use_underscore=1)
    ext_modules.append(ext)

    # c blas interface
    mod = 'cblas'
    sources = ['generic_blas1.pyf','generic_blas2.pyf','generic_blas3.pyf']
    ext = generic_extension(mod,sources,parent_package,use_underscore=0)
    ext_modules.append(ext)
    
    return ext_modules

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    #print configuration()
    setup(**configuration())
