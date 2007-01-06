#!/usr/bin/env python

from os.path import join
import sys

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse',parent_package,top_path)

    config.add_data_dir('tests')

    #config.add_library('sparsekit_src',
    #                   sources = [join('sparsekit','*.f')]
    #                   )

##    sources = ['spblas.f.src','spconv.f.src','sparsetools.pyf.src']
##    sources = [join('sparsetools',x) for x in sources]

##    config.add_extension('sparsetools',
##                         sources =  sources,
##                         )
    sources = ['sparsetools_wrap.cxx','sparsetools.py']
    sources = [join('sparsetools',x) for x in sources]

    config.add_extension('_sparsetools',
                         sources =  sources,
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
