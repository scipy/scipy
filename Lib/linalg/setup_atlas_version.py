#!/usr/bin/env python

import os
from distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict
from scipy_distutils.system_info import get_info,AtlasNotFoundError

def configuration (parent_package=''):
    package = 'linalg'
    config = default_config_dict(package,parent_package)
    del config['fortran_libraries']
    local_path = get_path(__name__)
    atlas_info = get_info('atlas_threads')
    if not atlas_info:
        atlas_info = get_info('atlas')
    if not atlas_info:
        raise AtlasNotFoundError,AtlasNotFoundError.__doc__
    ext = Extension('atlas_version',
                    sources=[os.path.join(local_path,'atlas_version.c')],
                    libraries=[atlas_info['libraries'][-1]],
                    library_dirs=atlas_info['library_dirs'])
    config['ext_modules'].append(ext)
    return config

if __name__ == '__main__':    
    from distutils.core import setup
    setup(**configuration())
