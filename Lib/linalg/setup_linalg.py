#!/usr/bin/env python
# To use:
#       python setup.py install
#

import os, sys, string
# this isn't in scipy_distutils -- need to figure out clean way to do that
from distutils import dep_util
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict
import interface_gen
from scipy_distutils.atlas_info import get_atlas_info

# needed now for pyf_extensions
import f2py2e

if os.name == 'nt':
    from scipy_distutils.mingw32_support import *

def configuration(parent_package=''):
    if parent_package:
        parent_package = parent_package + '.'    
    config = default_config_dict()
    config['packages'].append(parent_package+'linalg')
    config['packages'].append(parent_package+'linalg.tests')
    config['ext_modules'].extend(extensions(parent_package+'linalg'))
    return config
       
def generic_extension(mod_name,sources,parent_package='',use_underscore = 1):
    blas_libraries, lapack_libraries, atlas_library_dirs = get_atlas_info()
    
    if not atlas_library_dirs:
        msg = 'Atlas libraries not found.  Either install them in /usr/lib/atlas'\
              ' or /usr/local/lib/atlas and retry setup.py, or edit setup.py'\
              ' to specify your own blas and lapack directories and libs'
        raise ValueError, msg
    
    local_path = get_path(__name__)
    mod_file = mod_name + '.pyf'
    
    if parent_package and parent_package[-1] != '.':
        parent_package = parent_package + '.'    
    sources = [os.path.join(local_path,x) for x in sources]

    from distutils.dir_util import mkpath
    output_path = os.path.join('build','generated_pyfs')
    mkpath(output_path)
    mod_file = os.path.join(output_path,mod_file)

    # choose the
    gen_function = eval('interface_gen.generate_'+mod_name)
    
    if use_underscore:
        define_macros = []
    else:
        define_macros=[('NO_APPEND_FORTRAN',1)]
            
    if dep_util.newer_group(sources,mod_file):
        gen_function(local_path,output_path)
    
    print 'file:', mod_file    
    ext = Extension(parent_package+mod_name,[mod_file],
                     library_dirs=atlas_library_dirs,
                     libraries = lapack_libraries,
                     define_macros=define_macros,
                     f2py_options=['no-latex-doc',
                                   'no-makefile',
                                   'no-setup',
                                   'no-wrap-functions'])

    return ext                

def extensions(parent_package = ''):
    ext_modules = []
    
    # fortran lapack interface
    mod = 'flapack'
    sources = ['generic_lapack.pyf']
    ext = generic_extension(mod,sources,parent_package,use_underscore=1)
    ext_modules.append(ext)
    """
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
    """
    return ext_modules