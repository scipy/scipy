#!/usr/bin/env python

import os
from scipy.distutils.core import Extension
from scipy.distutils.misc_util import get_path, dot_join, Configuration

def configuration(parent_package='',parent_path=None):
    config = Configuration('io', parent_package, parent_path)
    config.add_extension('numpyio',
                         sources = ['numpyiomodule.c'])
    return config

if __name__ == '__main__':    
    from scipy.distutils.core import setup
    setup(**configuration(parent_path=''))
