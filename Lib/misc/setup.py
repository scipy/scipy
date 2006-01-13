
import os
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    config = Configuration('misc',parent_package, top_path)
    config.add_data_files('lena.dat')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration())
