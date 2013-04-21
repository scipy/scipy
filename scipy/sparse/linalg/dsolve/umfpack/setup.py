#!/usr/bin/env python
# 05.12.2005, c
# last change: 27.03.2006
from __future__ import division, print_function, absolute_import


def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, dict_append

    config = Configuration('umfpack', parent_package, top_path)
    config.add_data_dir('tests')

    umf_info = get_info('umfpack', notfound_action=1)

    umfpack_i_file = config.paths('umfpack.i')[0]

    def umfpack_i(ext, build_dir):
        if umf_info:
            return umfpack_i_file

    blas_info = get_info('blas_opt')
    build_info = {}
    dict_append(build_info, **umf_info)
    dict_append(build_info, **blas_info)

    config.add_extension('__umfpack',
                          sources=[umfpack_i],
                          depends=['umfpack.i'],
                          **build_info)

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
