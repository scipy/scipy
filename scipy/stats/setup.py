#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('stats', parent_package, top_path)

    config.add_data_dir('tests')

    statlib_src = [join('statlib', '*.f')]
    config.add_library('statlib', sources=statlib_src)

    # add statlib module
    config.add_extension('statlib',
        sources=['statlib.pyf'],
        f2py_options=['--no-wrap-functions'],
        libraries=['statlib'],
        depends=statlib_src
    )

    # add vonmises_cython module
    config.add_extension('vonmises_cython',
        sources=['vonmises_cython.c'],  # FIXME: use cython source,
        extra_link_args = ['-lpython2.7']
    )

    # add _rank module
    config.add_extension('_rank',
        sources=['_rank.c'],          # FIXME: use cython source
        extra_link_args = ['-lpython2.7']
    )

    # add futil module
    config.add_extension('futil',
        sources=['futil.f'],
    )

    # add mvn module
    config.add_extension('mvn',
        sources=['mvn.pyf','mvndst.f'],
    )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
