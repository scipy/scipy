#!/usr/bin/env python

import os
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path,Configuration,dot_join
join = os.path.join

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.system_info import get_info
    package = 'cluster'
    local_path = get_path(__name__,parent_path)
    config = Configuration(package,parent_package)

    config.add_extension('_vq',
        sources=[join('src', 'vq_wrap.cpp')])

    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(maintainer = "SciPy Developers",
          author = "Eric Jones",
          maintainer_email = "scipy-dev@scipy.org",
          description = "Clustering Algorithms (Information Theory)",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration()
          )
