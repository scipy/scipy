#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import os
from os.path import join


def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import get_sgemv_fix, get_g77_abi_wrappers, split_fortran_files

    config = Configuration('linalg',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[])
                      if k == 'ATLAS_INFO']+[None])[0]
    if atlas_version:
        print(('ATLAS version: %s' % atlas_version))

    # fblas:
    sources = ['fblas.pyf.src']
    sources += get_g77_abi_wrappers(lapack_opt)
    sources += get_sgemv_fix(lapack_opt)

    config.add_extension('_fblas',
                         sources=sources,
                         depends=['fblas_l?.pyf.src'],
                         extra_info=lapack_opt
                         )

    # flapack:
    sources = ['flapack.pyf.src']
    sources += get_g77_abi_wrappers(lapack_opt)

    config.add_extension('_flapack',
                         sources=sources,
                         depends=['flapack_user.pyf.src'],
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

    # _interpolative:
    routines_to_split = [
        'dfftb1',
        'dfftf1',
        'dffti1',
        'dsint1',
        'dzfft1',
        'id_srand',
        'idd_copyints',
        'idd_id2svd0',
        'idd_pairsamps',
        'idd_permute',
        'idd_permuter',
        'idd_random_transf0',
        'idd_random_transf0_inv',
        'idd_random_transf_init0',
        'idd_subselect',
        'iddp_asvd0',
        'iddp_rsvd0',
        'iddr_asvd0',
        'iddr_rsvd0',
        'idz_estrank0',
        'idz_id2svd0',
        'idz_permute',
        'idz_permuter',
        'idz_random_transf0_inv',
        'idz_random_transf_init0',
        'idz_random_transf_init00',
        'idz_realcomp',
        'idz_realcomplex',
        'idz_reco',
        'idz_subselect',
        'idzp_aid0',
        'idzp_aid1',
        'idzp_asvd0',
        'idzp_rsvd0',
        'idzr_asvd0',
        'idzr_reco',
        'idzr_rsvd0',
        'zfftb1',
        'zfftf1',
        'zffti1',
    ]
    print('Splitting linalg.interpolative Fortran source files')
    fnames = split_fortran_files(join(os.path.split(
                                          os.path.abspath(__file__))[0],
                                      'src', 'id_dist', 'src'),
                                 routines_to_split)
    fnames = [join('src', 'id_dist', 'src', f) for f in fnames]
    config.add_extension('_interpolative', fnames + ["interpolative.pyf"],
                         extra_info=lapack_opt
                         )

    # _calc_lwork:
    config.add_extension('_calc_lwork',
                         [join('src', 'calc_lwork.f')],
                         extra_info=lapack_opt)

    # _solve_toeplitz:
    config.add_extension('_solve_toeplitz',
                         sources=[('_solve_toeplitz.c')],
                         include_dirs=[get_numpy_include_dirs()])

    config.add_data_dir('tests')

    config.add_data_dir('benchmarks')

    config.add_data_files('blas_pointers.pxd')
    config.add_data_files('lapack_pointers.pxd')

    config.add_data_files('f2pyptr.h')

    sources = ['_blas_subroutine_wrappers.f']
    sources += get_g77_abi_wrappers(lapack_opt)
    sources += get_sgemv_fix(lapack_opt)
    config.add_extension('_blas_subroutine_wrappers',
                         sources=sources,
                         extra_info=lapack_opt)

    sources = ['_lapack_subroutine_wrappers.f']
    sources += get_g77_abi_wrappers(lapack_opt)
    config.add_extension('_lapack_subroutine_wrappers',
                         sources=sources,
                         extra_info=lapack_opt)

    config.add_extension('blas_pointers',
                         sources=['blas_pointers.c'],
                         include_dirs=['.'])

    config.add_extension('lapack_pointers',
                         sources=['lapack_pointers.c'],
                         include_dirs=['.'])

    config.add_extension('_cython_blas_wrappers',
                         sources=['_cython_blas_wrappers.c'],
                         include_dirs=['.'])

    config.add_extension('_cython_lapack_wrappers',
                         sources=['_cython_lapack_wrappers.c'],
                         include_dirs=['.'])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(top_path='').todict())
