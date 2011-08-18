#!/usr/bin/env python

import os
import sys
import re
from distutils.dep_util import newer_group, newer
from glob import glob
from os.path import join


def needs_veclib_wrapper(info):
    """Returns true if needs special veclib wrapper."""
    import re
    r_accel = re.compile("Accelerate")
    r_vec = re.compile("vecLib")
    res = False
    try:
        tmpstr = info['extra_link_args']
        for i in tmpstr:
            if r_accel.search(i) or r_vec.search(i):
                res = True
    except KeyError:
        pass

    return res


def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError

    from numpy.distutils.misc_util import Configuration

    config = Configuration('isolve',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    # iterative methods
    methods = ['BiCGREVCOM.f.src',
               'BiCGSTABREVCOM.f.src',
               'CGREVCOM.f.src',
               'CGSREVCOM.f.src',
#               'ChebyREVCOM.f.src',
               'GMRESREVCOM.f.src',
#               'JacobiREVCOM.f.src',
               'QMRREVCOM.f.src',
#               'SORREVCOM.f.src'
               ]

    if needs_veclib_wrapper(lapack_opt):
        methods += [join('FWRAPPERS', 'veclib_cabi_f.f'),
                    join('FWRAPPERS', 'veclib_cabi_c.c')]
    else:
        methods += [join('FWRAPPERS', 'dummy.f')]


    Util = ['STOPTEST2.f.src','getbreak.f.src']
    sources = Util + methods + ['_iterative.pyf.src']
    config.add_extension('_iterative',
                         sources=[join('iterative', x) for x in sources],
                         extra_info=lapack_opt,
                         depends=[join('iterative', 'FWRAPPERS', x) for x in
                         ['veclib_cabi_f.f', 'veclib_cabi_c.c', 'dummy.f']]
                         )

    config.add_data_dir('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
