#!/usr/bin/env python

import os
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path,default_config_dict,dot_join

def configuration(parent_package='',parent_path=None):
    package = 'cluster'
    local_path = get_path(__name__,parent_path)
    config = default_config_dict(package,parent_package)
    
    # This should really be fixed to use inline...
    sources = ['src/vq_wrap.cpp']
    sources = [os.path.join(local_path,x) for x in sources]

    ext = Extension(dot_join(parent_package,'cluster._vq'),sources)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(maintainer = "SciPy Developers",
          author = "eric jones",
          maintainer_email = "scipy-dev@scipy.org",
          description = "Clustering Algorithms (Information Theory)",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration()
          )
