#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('interpolate', parent_package, top_path)

    config.add_library('fitpack',
                       sources=[join('fitpack', '*.f')],
                      )

    config.add_extension('interpnd',
                         sources=['interpnd.c'])

    config.add_extension('_fitpack',
                         sources=['src/_fitpackmodule.c'],
                         libraries=['fitpack'],
                         depends = ['src/__fitpack.h','src/multipack.h']
                        )

    config.add_extension('dfitpack',
                         sources=['src/fitpack.pyf'],
                         libraries=['fitpack'],
                        )

    config.add_extension('_interpolate',
                         sources=['src/_interpolate.cpp'],
                         include_dirs = ['src'],
                         depends = ['src/interpolate.h'])

    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
