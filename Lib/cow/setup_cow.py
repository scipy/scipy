#!/usr/bin/env python

from scipy_distutils.misc_util import default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'cow'
    config = default_config_dict(package,parent_package)
    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
