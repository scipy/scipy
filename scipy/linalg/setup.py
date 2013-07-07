#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from os.path import join

from scipy._build_utils import needs_g77_abi_wrapper


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
        sources = ['fblas.pyf.src', join('src', 'fblaswrap.f')]

    # Note: `depends` needs to include fblaswrap(_veclib) for both files to be
    # included by "python setup.py sdist"
    config.add_extension('_fblas',
                         sources=sources,
                         depends=['fblas_l?.pyf.src',
                                    join('src', 'fblaswrap_veclib_c.c'),
                                    join('src', 'fblaswrap.f')],
                         extra_info=lapack_opt
                         )

    # flapack:
    config.add_extension('_flapack',
                         sources=['flapack.pyf.src'],
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

    # calc_lwork:
    config.add_extension('calc_lwork',
                         [join('src','calc_lwork.f')],
                         extra_info=lapack_opt
                         )

    # _interpolative:
    config.add_extension('_interpolative',
                         [join('src', 'id_dist', 'src', fn) for fn in 
                             ['dfft.f', 'idd_frm.f', 'idd_house.f', 'idd_id2svd.f', 'idd_id.f',
                                 'iddp_aid.f', 'iddp_asvd.f', 'iddp_rid.f', 'iddp_rsvd.f', 'idd_qrpiv.f',
                                 'iddr_aid.f', 'iddr_asvd.f', 'iddr_rid.f', 'iddr_rsvd.f', 'idd_sfft.f',
                                 'idd_snorm.f', 'idd_svd.f', 'id_rand.f', 'id_rtrans.f', 'idz_frm.f',
                                 'idz_house.f', 'idz_id2svd.f', 'idz_id.f', 'idzp_aid.f', 'idzp_asvd.f',
                                 'idzp_rid.f', 'idzp_rsvd.f', 'idz_qrpiv.f', 'idzr_aid.f', 'idzr_asvd.f',
                                 'idzr_rid.f', 'idzr_rsvd.f', 'idz_sfft.f', 'idz_snorm.f', 'idz_svd.f',
                                 ]
                             ] + ["interpolative.pyf"],
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
