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
    atlas_info = get_info('atlas')
    if atlas_info:
        ext = Extension('atlas_version',
                        sources=[os.path.join(local_path,'atlas_version.c')],
                        **atlas_info)
        config['ext_modules'].append(ext)
    else:
        print AtlasNotFoundError.__doc__
    return config

if __name__ == '__main__':    
    from distutils.core import setup
    setup(**configuration())
