#!/usr/bin/env python

import os
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'plt'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package,parent_package)

    data = ['lena.dat']
    data = [os.path.join(local_path,x) for x in data]               
    config['data_files'].extend( [(os.path.join(parent_package,'plt'),data)])
  
    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
