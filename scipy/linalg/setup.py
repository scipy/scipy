#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
from os.path import join

from scipy._build_utils import needs_g77_abi_wrapper, split_fortran_files


def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError

    from numpy.distutils.misc_util import Configuration

    config = Configuration('linalg',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[])
                      if k == 'ATLAS_INFO']+[None])[0]
    if atlas_version:
        print(('ATLAS version: %s' % atlas_version))

    # fblas:
    if needs_g77_abi_wrapper(lapack_opt):
        sources = ['fblas.pyf.src', join('src', 'fblaswrap_veclib_c.c')],
    else:
        sources = ['fblas.pyf.src', join('src', 'fblaswrap_dummy.f')]

    # Note: `depends` needs to include fblaswrap(_veclib) for both files to be
    # included by "python setup.py sdist"
    config.add_extension('_fblas',
                         sources=sources,
                         depends=['fblas_l?.pyf.src',
                                  join('src', 'fblaswrap_veclib_c.c'),
                                  join('src', 'fblaswrap_dummy.f')],
                         extra_info=lapack_opt
                         )

    # flapack:
    if needs_g77_abi_wrapper(lapack_opt):
        sources = ['flapack.pyf.src', join('src', 'flapackwrap_veclib.f')],
    else:
        sources = ['flapack.pyf.src', join('src', 'flapackwrap_dummy.f')]

    config.add_extension('_flapack',
                         sources=sources,
                         depends=['flapack_user.pyf.src',
                                  join('src', 'flapackwrap_veclib.f'),
                                  join('src', 'flapackwrap_dummy.f')],
                         extra_info=lapack_opt
                         )

    if atlas_version is not None:
        # cblas:
        config.add_extension('_cblas',
                             sources=['cblas.pyf.src'],
                             depends=['cblas.pyf.src', 'cblas_l1.pyf.src'],
                             extra_info=lapack_opt
                             )

        # clapack:
        config.add_extension('_clapack',
                             sources=['clapack.pyf.src'],
                             depends=['clapack.pyf.src'],
                             extra_info=lapack_opt
                             )

    # _flinalg:
    config.add_extension('_flinalg',
                         sources=[join('src','det.f'),join('src','lu.f')],
                         extra_info=lapack_opt
                         )

    # calc_lwork:
    config.add_extension('calc_lwork',
                         [join('src','calc_lwork.f')],
                         extra_info=lapack_opt
                         )

    # _interpolative:
    print('Splitting linalg.interpolative Fortran source files')
    fnames = split_fortran_files(join(os.path.split(
                                          os.path.abspath(__file__))[0],
                                      'src', 'id_dist', 'src'))
    fnames = [join('src', 'id_dist', 'src', f) for f in fnames]
    config.add_extension('_interpolative', fnames + ["interpolative.pyf"],
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
