#!/usr/bin/env python

import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, dot_join, default_config_dict

def configuration(parent_package='',parent_path=None):
    from scipy_distutils.system_info import get_info, dict_append
    package = 'io'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package,parent_package)

    numpy_info = get_info('numpy',notfound_action=2)

    sources = ['numpyiomodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
    ext = Extension(dot_join(parent_package,package,'numpyio'),sources,
                    **numpy_info)
    config['ext_modules'].append(ext)
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
