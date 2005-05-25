#!/usr/bin/env python

from glob import glob
import os

def configuration(parent_package='', parent_path=None):
    from scipy_distutils.system_info import get_info, dict_append
    from scipy_distutils.misc_util import default_config_dict, \
         dot_join, get_path    

    package = 'image'
    config = default_config_dict(package,parent_package)

    local_path = get_path(__name__, parent_path)
    image_path = os.path.join(parent_package,'image')

    color_files = glob(os.path.join(local_path, '*.txt'))
    data_path = os.path.join(image_path, 'colordata')
    config['data_files'].extend([(data_path, color_files)])

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
