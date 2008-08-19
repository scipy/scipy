"""setup"""

import os, sys, setuptools

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)
    config.add_extension('_triang',
                                    sources = ['_triang.cpp'],
                                    depends = [],
                                    )
    
    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        **configuration().todict()
    )
