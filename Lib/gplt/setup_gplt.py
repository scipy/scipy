#!/usr/bin/env python

import os, sys
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'gplt'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package,parent_package)

    if sys.platform == 'win32':
        data = ['wgnuplot.exe', 'gnuplot_helper.exe']
        data = [os.path.join(local_path,x) for x in data]
        config['data_files'].append((os.path.join(parent_package,'gplt'),data))
    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
