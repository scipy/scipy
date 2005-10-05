#!/usr/bin/env python

import os
from scipy.distutils.misc_util import get_path, default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'ga'
    config = default_config_dict(package,parent_package)
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(**configuration())
