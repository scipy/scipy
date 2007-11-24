#!/usr/bin/env python

import os
import sys
from os.path import join
from distutils.sysconfig import get_python_inc

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('special', parent_package, top_path)

    config.add_sconscript('SConstruct')
    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
