from __future__ import division, print_function, absolute_import

import os
import re
import glob
from os.path import join


def patch_int64_pyf(filename, basename_suffix='_64'):
    """
    Patch .pyf files to ILP64 ABI.

    Does the replacements::

    - ``npy_int -> npy_int64``
    - Appends `basename_suffix` to module name and file basename.
    - Appends `basename_suffix` to include directives.

    """
    from scipy._build_utils import write_file_content

    with open(filename, 'r') as f:
        text = f.read()

    text = re.sub(r"npy_int(?!\w)", "npy_int64", text, flags=re.I)
    text = re.sub(r"python module ([a-z0-9_]+)",
                  r"python module \1{}".format(basename_suffix),
                  text, flags=re.I)
    text = re.sub(r"include\s+'([a-z0-9_]+)(\.[a-z0-9.]+)'",
                  r"include '\1{}\2'".format(basename_suffix),
                  text, flags=re.I)

    basename = os.path.basename(filename)
    if basename.endswith('.pyf.src'):
        base, ext = basename[:-8], basename[-8:]
    else:
        base, ext = os.path.splitext(basename)

    ext_dir = os.path.relpath(os.path.dirname(filename), os.getcwd())
    dst = os.path.join('build', 'build_src.common', ext_dir,
                       base + basename_suffix + ext)

    write_file_content(dst, text)

    return os.path.relpath(dst, os.path.dirname(filename))


def patch_int64_pyf_glob(filenames, basename_suffix='_64'):
    new_filenames = []
    for fn_glob in filenames:
        fn_glob = os.path.join(os.path.dirname(__file__), fn_glob)
        for fn in glob.glob(fn_glob):
            if not (fn.endswith('.pyf') or fn.endswith('.pyf.src')):
                new_filenames.append(fn)
                continue

            new_filenames.append(patch_int64_pyf(fn, basename_suffix))

    return new_filenames


def configuration(parent_package='', top_path=None):
    from distutils.sysconfig import get_python_inc
    from scipy._build_utils.system_info import get_info, NotFoundError, numpy_info
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import (get_g77_abi_wrappers, uses_blas64,
                                    blas_ilp64_pre_build_hook, get_f2py_int64_options)

    config = Configuration('linalg', parent_package, top_path)

    lapack_opt = get_info('lapack_opt')

    atlas_version = ([v[3:-3] for k, v in lapack_opt.get('define_macros', [])
                      if k == 'ATLAS_INFO']+[None])[0]
    if atlas_version:
        print(('ATLAS version: %s' % atlas_version))

    if uses_blas64():
        lapack_ilp64_opt = get_info('lapack_ilp64_opt', 2)

    # fblas:
    sources = ['fblas.pyf.src']
    sources += get_g77_abi_wrappers(lapack_opt)

    config.add_extension('_fblas',
                         sources=sources,
                         depends=['fblas_l?.pyf.src'],
                         extra_info=lapack_opt
                         )

    if uses_blas64():
        ext = config.add_extension('_fblas_64',
                                   sources=patch_int64_pyf_glob(sources),
                                   depends=patch_int64_pyf_glob(['fblas_l?.pyf.src']),
                                   extra_info=lapack_ilp64_opt,
                                   f2py_options=get_f2py_int64_options())
        ext._pre_build_hook = blas_ilp64_pre_build_hook(lapack_ilp64_opt)

    # flapack:
    sources = ['flapack.pyf.src']
    sources += get_g77_abi_wrappers(lapack_opt)
    dep_pfx = join('src', 'lapack_deprecations')
    deprecated_lapack_routines = [join(dep_pfx, c + 'gegv.f') for c in 'cdsz']
    sources += deprecated_lapack_routines

    config.add_extension('_flapack',
                         sources=sources,
                         depends=['flapack_gen.pyf.src',
                                  'flapack_gen_banded.pyf.src',
                                  'flapack_gen_tri.pyf.src',
                                  'flapack_pos_def.pyf.src',
                                  'flapack_pos_def_tri.pyf.src',
                                  'flapack_sym_herm.pyf.src',
                                  'flapack_other.pyf.src',
                                  'flapack_user.pyf.src'],
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
                         sources=[join('src', 'det.f'), join('src', 'lu.f')],
                         extra_info=lapack_opt
                         )

    # _interpolative:
    config.add_extension('_interpolative',
                         sources=[join('src', 'id_dist', 'src', '*.f'),
                                  "interpolative.pyf"],
                         extra_info=lapack_opt
                         )

    # _solve_toeplitz:
    config.add_extension('_solve_toeplitz',
                         sources=[('_solve_toeplitz.c')],
                         include_dirs=[get_numpy_include_dirs()])

    config.add_data_dir('tests')

    # Cython BLAS/LAPACK
    config.add_data_files('cython_blas.pxd')
    config.add_data_files('cython_lapack.pxd')

    sources = ['_blas_subroutine_wrappers.f', '_lapack_subroutine_wrappers.f']
    sources += get_g77_abi_wrappers(lapack_opt)
    includes = numpy_info().get_include_dirs() + [get_python_inc()]
    config.add_library('fwrappers', sources=sources, include_dirs=includes)

    config.add_extension('cython_blas',
                         sources=['cython_blas.c'],
                         depends=['cython_blas.pyx', 'cython_blas.pxd',
                                  'fortran_defs.h', '_blas_subroutines.h'],
                         include_dirs=['.'],
                         libraries=['fwrappers'],
                         extra_info=lapack_opt)

    config.add_extension('cython_lapack',
                         sources=['cython_lapack.c'],
                         depends=['cython_lapack.pyx', 'cython_lapack.pxd',
                                  'fortran_defs.h', '_lapack_subroutines.h'],
                         include_dirs=['.'],
                         libraries=['fwrappers'],
                         extra_info=lapack_opt)

    config.add_extension('_decomp_update',
                         sources=['_decomp_update.c'])

    # Add any license files
    config.add_data_files('src/id_dist/doc/doc.tex')
    config.add_data_files('src/lapack_deprecations/LICENSE')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    from linalg_version import linalg_version

    setup(version=linalg_version,
          **configuration(top_path='').todict())
