#!/usr/bin/env python

from os.path import join
import sys

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparse',parent_package,top_path)

    config.add_data_dir('tests')
    config.add_subpackage('linalg')

    # Adding a Python file as a "source" file for an extension is something of
    # a hack, but it works to put it in the right place.
    sources = [join('sparsetools', x) for x in
        ['sparsetools.py', 'sparsetools_wrap.cxx']]
    config.add_extension('_sparsetools',
                         sources=sources,
                         include_dirs=['sparsetools'],
                         )
##     sparsetools_i_file = config.paths(join('sparsetools','sparsetools.i'))[0]
##     def sparsetools_i(ext, build_dir):
##             return sparsetools_i_file
##     config.add_extension('_sparsetools',
##                          sources= [sparsetools_i_file],
##                          include_dirs=['sparsetools'],
##                          depends = [join('sparsetools', x) for x in
##                                     ['sparsetools.i', 'sparsetools.h']]
##                          )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
