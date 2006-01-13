#!/usr/bin/env python

import os
from numpy.distutils.misc_util import get_path, default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'plt'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package,parent_package)
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration())
