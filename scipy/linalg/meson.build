__init__py = fs.copyfile('__init__.py')
_cy_array_utils_pxd = fs.copyfile('_cythonized_array_utils.pxd')

_cython_tree += fs.copyfile('__init__.pxd')

_generate_pyx = find_program('_generate_pyx.py')

ilp64_flag = cython_blas_ilp64 ? '--ilp64' : []
cython_linalg = custom_target('cython_linalg',
  output: [
    'cython_blas.pyx',
    'cython_lapack.pyx',
    'cython_blas.pxd',
    'cython_lapack.pxd',
    '_blas_subroutines.h',
    '_lapack_subroutines.h'
  ],
  input: '_generate_pyx.py',
  command: [_generate_pyx, '-o', '@OUTDIR@', accelerate_flag, ilp64_flag],
  depend_files: [
    'cython_blas_signatures.txt',
    'cython_lapack_signatures.txt',
    '../_build_utils/_generate_blas_wrapper.py',
    '../_build_utils/_wrappers_common.py',
  ],
  install: true,
  # Install only the .pxd files
  install_dir: [
    false,
    false,
    py3.get_install_dir() / 'scipy/linalg',
    py3.get_install_dir() / 'scipy/linalg',
    false,
    false,
  ],
  install_tag: 'devel'
)
cython_blas_pxd = cython_linalg[2]
cython_lapack_pxd = cython_linalg[3]

# pyx -> c, pyx -> cpp generators, depending on __init__.py here.
linalg_init_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, __init__py])

# pyx -> c, pyx -> cpp generators, depending on _cythonized_array_utils.pxd
linalg_init_utils_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, __init__py, _cy_array_utils_pxd, cython_lapack_pxd])

# pyx -> c, pyx -> cpp generators, depending on copied pxd files and init
linalg_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, __init__py, cython_blas_pxd])#, cython_linalg])

link_language_blas_lapack = 'fortran'
if get_option('_without-fortran')
  link_language_blas_lapack = 'c'
endif

if needs_lp64_fblas
  fblas_module = custom_target('fblas_module',
    output: ['_fblasmodule.c'],
    input: 'fblas.pyf.src',
    command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@'] + f2py_freethreading_arg,
    depend_files:
      [
        'fblas_l1.pyf.src',
        'fblas_l2.pyf.src',
        'fblas_l3.pyf.src',
      ]
  )

  # Note: we're linking LAPACK on purpose here. For some routines (e.g., spmv)
  # the float/double routines are in BLAS while the complex routines are in
  # LAPACK - we have historically put these in `_fblas`.
  py3.extension_module('_fblas',
    fblas_module,
    link_args: version_link_args,
    dependencies: [lapack_lp64_dep, fortranobject_dep],
    link_with: blas_lapack_wrapper_lib_lp64,
    install: true,
    link_language: link_language_blas_lapack,
    subdir: 'scipy/linalg'
  )

  flapack_module = custom_target('flapack_module',
    output: ['_flapackmodule.c'],
    input: 'flapack.pyf.src',
    command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@'] + f2py_freethreading_arg,
    depend_files:
      [
        'flapack_user.pyf.src',
        'flapack_gen.pyf.src',
        'flapack_gen_banded.pyf.src',
        'flapack_gen_tri.pyf.src',
        'flapack_sym_herm.pyf.src',
        'flapack_pos_def.pyf.src',
        'flapack_pos_def_tri.pyf.src',
        'flapack_other.pyf.src',
      ]
  )

  # Note that -Wno-empty-body is Clang-specific and comes from `callstatement`s
  # in flapack_other.pyf.src
  py3.extension_module('_flapack',
    flapack_module,
    c_args: [Wno_empty_body],
    link_args: version_link_args,
    dependencies: [lapack_lp64_dep, fortranobject_dep],
    link_with: blas_lapack_wrapper_lib_lp64,
    install: true,
    subdir: 'scipy/linalg'
  )
endif

if use_ilp64
  fblas64_module = custom_target('fblas64_module',
    output: ['_fblas_64module.c'],
    input: 'fblas_64.pyf.src',
    command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@'] + f2py_ilp64_opts + f2py_freethreading_arg,
    depend_files:
      [
        'fblas_l1.pyf.src',
        'fblas_l2.pyf.src',
        'fblas_l3.pyf.src',
        int64_f2cmap,
      ]
  )

  py3.extension_module('_fblas_64',
    fblas64_module,
    link_args: version_link_args,
    dependencies: [lapack_dep, fortranobject_dep],
    link_with: blas_lapack_wrapper_lib_ilp64,
    install: true,
    link_language: link_language_blas_lapack,
    subdir: 'scipy/linalg'
  )

  flapack64_module = custom_target('flapack64_module',
    output: ['_flapack_64module.c'],
    input: 'flapack_64.pyf.src',
    command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@'] + f2py_ilp64_opts + f2py_freethreading_arg,
    depend_files:
      [
        'flapack_user.pyf.src',
        'flapack_gen.pyf.src',
        'flapack_gen_banded.pyf.src',
        'flapack_gen_tri.pyf.src',
        'flapack_sym_herm.pyf.src',
        'flapack_pos_def.pyf.src',
        'flapack_pos_def_tri.pyf.src',
        'flapack_other.pyf.src',
        int64_f2cmap,
      ]
  )

  py3.extension_module('_flapack_64',
    flapack64_module,
    c_args: [Wno_empty_body],
    link_args: version_link_args,
    dependencies: [lapack_dep, fortranobject_dep],
    link_with: blas_lapack_wrapper_lib_ilp64,
    install: true,
    link_language: link_language_blas_lapack,
    subdir: 'scipy/linalg'
  )
endif

# _decomp_interpolative
py3.extension_module('_decomp_interpolative',
  linalg_cython_gen.process('_decomp_interpolative.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

# batched linalg inv module
py3.extension_module('_batched_linalg',
  [
    'src/_common_array_utils.hh',
    'src/_linalg_lu_det.hh',
    'src/_linalg_inv.hh',
    'src/_linalg_solve.hh',
    'src/_npymath.hh',
    'src/_batched_linalg_module.cc',
  ],
  include_directories: ['src'],
  dependencies: [np_dep, lapack_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

# _solve_toeplitz
py3.extension_module('_solve_toeplitz',
  linalg_init_cython_gen.process('_solve_toeplitz.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

# _matfuncs_sqrtm_triu:
py3.extension_module('_matfuncs_sqrtm_triu',
  linalg_init_cython_gen.process('_matfuncs_sqrtm_triu.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

cython_blas_lapack_dep = cython_blas_ilp64 ? lapack_dep : lapack_lp64_dep
cython_blas = py3.extension_module('cython_blas',
  [
    linalg_cython_gen.process(cython_linalg[0]),  # cython_blas.pyx
    cython_linalg[4],  # _blas_subroutines.h
  ],
  c_args: cython_c_args,
  link_args: version_link_args,
  dependencies: [cython_blas_lapack_dep, np_dep],
  install: true,
  subdir: 'scipy/linalg'
)

cython_lapack = py3.extension_module('cython_lapack',
  [
    linalg_cython_gen.process(cython_linalg[1]),  # cython_lapack.pyx
    cython_linalg[5],  # _lapack_subroutines.h
  ],
  c_args: cython_c_args,
  link_args: version_link_args,
  dependencies: [cython_blas_lapack_dep, np_dep],
  install: true,
  subdir: 'scipy/linalg'
)

_decomp_update_pyx = custom_target('_decomp_update',
  output: '_decomp_update.pyx',
  input: '_decomp_update.pyx.in',
  command: [tempita, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_decomp_update',
  linalg_cython_gen.process(_decomp_update_pyx),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

py3.extension_module('_internal_matfuncs',
  [
    '_matfuncsmodule.c',
    'src/_matfuncs_expm.c',
    'src/_matfuncs_sqrtm.c',
  ],
  include_directories: ['src/'],
  dependencies: [np_dep, lapack_dep],
  install: true,
  subdir: 'scipy/linalg'
)


_cythonized_array_utils = py3.extension_module('_cythonized_array_utils',
  linalg_init_utils_cython_gen.process('_cythonized_array_utils.pyx'),
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/linalg'
)

if use_pythran
  py3.extension_module('_linalg_pythran',
    pythran_gen.process('_linalg_pythran.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/linalg'
  )
else
  py3.install_sources(
    ['_linalg_pythran.py'],
    subdir: 'scipy/linalg'
  )
endif


python_sources = [
  '__init__.py',
  '__init__.pxd',
  '_basic.py',
  '_cythonized_array_utils.pxd',
  '_cythonized_array_utils.pyi',
  '_decomp.py',
  '_decomp_cholesky.py',
  '_decomp_cossin.py',
  '_decomp_ldl.py',
  '_decomp_lu.py',
  '_decomp_polar.py',
  '_decomp_qr.py',
  '_decomp_qz.py',
  '_decomp_schur.py',
  '_decomp_svd.py',
  '_expm_frechet.py',
  '_matfuncs.py',
  '_matfuncs_expm.pyi',
  '_matfuncs_inv_ssq.py',
  '_matfuncs_sqrtm.py',
  '_misc.py',
  '_procrustes.py',
  '_sketches.py',
  '_solvers.py',
  '_special_matrices.py',
  '_testutils.py',
  'basic.py',
  'blas.py',
  'decomp.py',
  'decomp_cholesky.py',
  'decomp_lu.py',
  'decomp_qr.py',
  'decomp_schur.py',
  'decomp_svd.py',
  'interpolative.py',
  'lapack.py',
  'matfuncs.py',
  'misc.py',
  'special_matrices.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/linalg'
)

subdir('tests')
