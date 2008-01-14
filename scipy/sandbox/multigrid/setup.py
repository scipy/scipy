#!/usr/bin/env python

from os.path import join
import sys

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('multigrid',parent_package,top_path)

    config.add_subpackage('gallery')

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','sample_data'))

    # Adding a Python file as a "source" file for an extension is something of
    # a hack, but it works to put it in the right place.
    sources = [join('multigridtools', x) for x in ['multigridtools.py', 'multigridtools_wrap.cxx']]
    config.add_extension('_multigridtools',
                         sources=sources,
                         include_dirs=['multigridtools'])
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
