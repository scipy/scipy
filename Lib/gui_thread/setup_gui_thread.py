#!/usr/bin/env python

import os
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package='',parent_path=None):
    package = 'gui_thread'
    config = default_config_dict(package,parent_package)
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(maintainer = "SciPy Developers",
          maintainer_email = "scipy-dev@scipy.org",
          description = "SciPy test module",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration())
