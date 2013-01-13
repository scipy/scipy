#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
from distutils.dep_util import newer_group, newer
from os.path import join

from scipy._build_utils import needs_g77_abi_wrapper

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError

    from numpy.distutils.misc_util import Configuration

    config = Configuration('linalg',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print(('ATLAS version: %s' % atlas_version))

    target_dir = ''

    # fblas:
    if needs_g77_abi_wrapper(lapack_opt):
        sources = ['fblas.pyf.src', join('src', 'fblaswrap_veclib_c.c')],
    else:
        sources = ['fblas.pyf.src', join('src', 'fblaswrap.f')]

    # Note: `depends` needs to include fblaswrap(_veclib) for both files to be
    # included by "python setup.py sdist"
    config.add_extension('_fblas',
                         sources = sources,
                         depends = ['fblas_l?.pyf.src',
                                    join('src', 'fblaswrap_veclib_c.c'),
                                    join('src', 'fblaswrap.f')],
                         extra_info = lapack_opt
                         )

    # flapack:
    config.add_extension('_flapack',
                         sources = ['flapack.pyf.src'],
                         depends = ['flapack_user.pyf.src'],
                         extra_info = lapack_opt
                         )

    if atlas_version is not None:
        # cblas:
        config.add_extension('_cblas',
                             sources = ['cblas.pyf.src'],
                             depends = ['cblas.pyf.src', 'cblas_l1.pyf.src'],
                             extra_info = lapack_opt
                             )

        # clapack:
        config.add_extension('_clapack',
                             sources = ['clapack.pyf.src'],
                             depends = ['clapack.pyf.src'],
                             extra_info = lapack_opt
                             )

    # _flinalg:
    config.add_extension('_flinalg',
                         sources = [join('src','det.f'),join('src','lu.f')],
                         extra_info = lapack_opt
                         )

    # calc_lwork:
    config.add_extension('calc_lwork',
                         [join('src','calc_lwork.f')],
                         extra_info = lapack_opt
                         )

    config.add_data_dir('tests')
    config.add_data_dir('benchmarks')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(top_path='').todict())
