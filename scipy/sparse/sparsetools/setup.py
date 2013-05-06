#!/usr/bin/env python
from __future__ import division, print_function, absolute_import


def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('sparsetools',parent_package,top_path)

    for fmt in ['csr','csc','coo','bsr','dia','csgraph']:
        sources = [fmt + '_wrap.cxx']
        depends = [fmt + '.h']
        config.add_extension('_' + fmt, sources=sources,
            define_macros=[('__STDC_FORMAT_MACROS', 1)],
            depends=depends)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
