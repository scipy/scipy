from __future__ import division, print_function, absolute_import

import os
import re
import glob
from os.path import join


def patch_int64_pyf(filename, lapack_opt, basename_suffix='_64'):
    """
    Patch .pyf files to ILP64 ABI.

    Does the replacements::

    - ``int -> npy_int64``
    - Appends `basename_suffix` to file basename, module names,
      and callback identifier names.

    """
    from scipy._build_utils import write_file_content

    with open(filename, 'r') as f:
        text = f.read()

    text = re.sub(r"^int(?!\w)", r"npy_int64", text, flags=re.I)
    text = re.sub(r"([^\w])int(?!\w)", r"\1npy_int64", text, flags=re.I)
    text = re.sub(r"python\s+module\s+([a-z0-9_]+)",
                  r"python module \1{}".format(basename_suffix),
                  text, flags=re.I)
    text = re.sub(r"([^\w])use ([a-z0-9_]+)",
                  r"\1use \2{}".format(basename_suffix),
                  text, flags=re.I)
    text = re.sub(r"([^\w])(cb_[a-z_0-9<>]+__user__routines)(?!\w)",
                  r"\1\2{}".format(basename_suffix),
                  text, flags=re.I)
    text = re.sub(r"include\s+'([a-z0-9_]+)(\.[a-z0-9.]+)'",
                  r"include '\1{}\2'".format(basename_suffix),
                  text, flags=re.I)

    macros = dict(lapack_opt.get('define_macros', []))
    suffix = macros.get('BLAS_SYMBOL_SUFFIX', '')
    prefix = macros.get('BLAS_SYMBOL_PREFIX', '')
    if suffix or prefix:
        # Add rename macro declarations in C code
        text = re.sub(r"^(\s*python module.*)$",
                      ("\\1\n"
                       "usercode '''\n"
                       "  #include \"blas64-prefix-defines.h\"\n"
                       "'''\n\n"),
                      text, flags=re.I | re.M)

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


def patch_int64_pyf_glob(filenames, lapack_opt, basename_suffix='_64'):
    new_filenames = []
    for fn_glob in filenames:
        fn_glob = os.path.join(os.path.dirname(__file__), fn_glob)
        for fn in glob.glob(fn_glob):
            if not (fn.endswith('.pyf') or fn.endswith('.pyf.src')):
                new_filenames.append(fn)
                continue

            new_filenames.append(patch_int64_pyf(fn, lapack_opt, basename_suffix))

    return new_filenames


def configuration(parent_package='', top_path=None):
    from distutils.sysconfig import get_python_inc
    from scipy._build_utils.system_info import get_info, NotFoundError, numpy_info
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import (get_g77_abi_wrappers, uses_blas64,
                                    blas_ilp64_pre_build_hook, get_f2py_int64_options,
                                    generic_pre_build_hook)

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
                                   sources=patch_int64_pyf_glob(sources, lapack_ilp64_opt),
                                   depends=patch_int64_pyf_glob(['fblas_l?.pyf.src'], lapack_ilp64_opt),
                                   extra_info=lapack_ilp64_opt,
                                   f2py_options=get_f2py_int64_options())
        ext._pre_build_hook = blas_ilp64_pre_build_hook(lapack_ilp64_opt)

    # flapack:
    sources = ['flapack.pyf.src']
    sources += get_g77_abi_wrappers(lapack_opt)
    dep_pfx = join('src', 'lapack_deprecations')
    deprecated_lapack_routines = [join(dep_pfx, c + 'gegv.f') for c in 'cdsz']
    sources += deprecated_lapack_routines
    flapack_depends = ['flapack_gen.pyf.src',
                       'flapack_gen_banded.pyf.src',
                       'flapack_gen_tri.pyf.src',
                       'flapack_pos_def.pyf.src',
                       'flapack_pos_def_tri.pyf.src',
                       'flapack_sym_herm.pyf.src',
                       'flapack_other.pyf.src',
                       'flapack_user.pyf.src']

    config.add_extension('_flapack',
                         sources=sources,
                         depends=flapack_depends,
                         extra_info=lapack_opt
                         )

    if uses_blas64():
        ext = config.add_extension('_flapack_64',
                                   sources=patch_int64_pyf_glob(sources, lapack_ilp64_opt),
                                   depends=patch_int64_pyf_glob(flapack_depends, lapack_ilp64_opt),
                                   extra_info=lapack_ilp64_opt,
                                   f2py_options=get_f2py_int64_options())
        ext._pre_build_hook = blas_ilp64_pre_build_hook(lapack_ilp64_opt)

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
    id_dist_src = [join('src', 'id_dist', 'src', '*.f')]

    if uses_blas64():
        # Patch id_dist sources on the fly, changing *gesdd -> *gesdd_id_dist,
        # and include LAPACK wrappers that define it.
        def id_dist_pre_build_hook(cmd, ext):
            fcompiler_flags = {'gnu95': ['-cpp', '-ffree-line-length-none', '-ffixed-line-length-none']}

            macros = dict(lapack_ilp64_opt.get('define_macros', []))
            prefix = macros.get('BLAS_SYMBOL_PREFIX', '')
            suffix = macros.get('BLAS_SYMBOL_SUFFIX', '')
            if suffix:
                if not suffix.endswith('_'):
                    raise ValueError("Incompatible BLAS symbol suffix (has to end in _)")
                suffix = '_' + suffix[:-1]

            def patch_source(filename, old_text):
                if os.path.basename(filename) == 'id_dist_gesdd.f':
                    text = ("#define dgesdd {0}dgesdd{1}\n"
                            "#define zgesdd {0}zgesdd{1}\n").format(prefix, suffix)
                else:
                    text = ("#define dgesdd dgesdd_id_dist\n"
                            "#define zgesdd zgesdd_id_dist\n")
                text += old_text
                return text

            return generic_pre_build_hook(cmd, ext,
                                          fcompiler_flags=fcompiler_flags,
                                          patch_source_func=patch_source,
                                          source_fnpart='_id_dist_64')

        id_dist_src.append(os.path.join('src', 'id_dist_gesdd.f'))
        id_dist_lapack_opt = lapack_ilp64_opt
    else:
        id_dist_pre_build_hook = None
        id_dist_lapack_opt = lapack_opt

    config.add_library('id_dist_src',
                       sources=id_dist_src,
                       _pre_build_hook=id_dist_pre_build_hook)

    ext = config.add_extension('_interpolative',
                               sources=["interpolative.pyf"],
                               libraries=['id_dist_src'],
                               extra_info=id_dist_lapack_opt)
    ext._pre_build_hook = id_dist_pre_build_hook

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
