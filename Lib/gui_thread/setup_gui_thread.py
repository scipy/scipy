#!/usr/bin/env python

import os
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join

def configuration(parent_package=''):
    config = default_config_dict('gui_thread',parent_package)
    local_path = get_path(__name__)
    test_path = os.path.join(local_path,'tests')
    config['packages'].append(dot_join(parent_package,'gui_thread.tests')) 
    config['package_dir']['gui_thread.tests'] = test_path
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
