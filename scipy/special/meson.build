# TODO: 64-bit BLAS and LAPACK

_ufuncs_pxi_pxd_sources = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_agm.pxd'),
  fs.copyfile('_boxcox.pxd'),
  fs.copyfile('_cdflib_wrappers.pxd'),
  fs.copyfile('_complexstuff.pxd'),
  fs.copyfile('_convex_analysis.pxd'),
  fs.copyfile('_cunity.pxd'),
  fs.copyfile('_ellip_harm.pxd'),
  fs.copyfile('_ellip_harm_2.pxd'),
  fs.copyfile('_ellipk.pxd'),
  fs.copyfile('_factorial.pxd'),
  fs.copyfile('_hyp0f1.pxd'),
  fs.copyfile('_hypergeometric.pxd'),
  fs.copyfile('_legacy.pxd'),
  fs.copyfile('_ndtri_exp.pxd'),
  fs.copyfile('_sici.pxd'),
  fs.copyfile('_spence.pxd'),
  fs.copyfile('_xlogy.pxd'),
  fs.copyfile('orthogonal_eval.pxd'),
  fs.copyfile('sf_error.pxd'),
  fs.copyfile('_ufuncs_extra_code.pxi'),
  fs.copyfile('_ufuncs_extra_code_common.pxi'),
]

if meson.get_compiler('cpp').get_id() in ['msvc', 'clang-cl', 'intel-cl']
  scipy_export_dll_args = ['-DSCIPY_DLL_EXPORTS']
  scipy_import_dll_args = ['-DSCIPY_DLL_IMPORTS']
else
  scipy_export_dll_args = []
  scipy_import_dll_args = []
endif

sf_error_state_lib = shared_library('sf_error_state',
  ['sf_error_state.c'],
  include_directories: ['../_lib', '../_build_utils/src'],
  c_args: scipy_export_dll_args,
  install: true,
  install_dir: py3.get_install_dir() / 'scipy/special',
)

sf_error_state_dep = declare_dependency(
  compile_args: ['-DSP_SPECFUN_ERROR'] + scipy_import_dll_args,
  link_with: sf_error_state_lib,
)

ufuncs_sources = [
  '_cosine.c',
  'special_wrappers.cpp',
  'sf_error.cc',
  'dd_real_wrappers.cpp'
]

ufuncs_cxx_sources = [
  '_faddeeva.cxx',
  '_wright.cxx',
  'ellint_carlson_wrap.cxx',
  'Faddeeva.cc',
  'sf_error.cc',
  'wright.cc'
]

cdflib_lib = static_library('cdflib',
  'cdflib.c',
  include_directories: ['../_lib', '../_build_utils/src'],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('_special_ufuncs',
  ['_special_ufuncs.cpp', '_special_ufuncs_docs.cpp', 'sf_error.cc'],
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [np_dep, sf_error_state_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN'
)

py3.extension_module('_gufuncs',
  ['_gufuncs.cpp', '_gufuncs_docs.cpp', 'sf_error.cc'],
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [np_dep, sf_error_state_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN'
)

py3.extension_module('_specfun',
  [cython_gen_cpp.process('_specfun.pyx')],
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special'
)

# TODO: this installs all generated files, while we want to install only:
#       _ufuncs.pyi, cython_special.pxd
#       There's no good way around this currently, because py3.install_sources
#       doesn't accept generated targets. See TODO near the end of
#       scipy/linalg/meson.build for more details.
_generate_pyx = find_program('_generate_pyx.py')
cython_special = custom_target('cython_special',
  output: [
    '_ufuncs.pyx',
    '_ufuncs_defs.h',
    '_ufuncs_cxx.pyx',
    '_ufuncs_cxx.pxd',
    '_ufuncs_cxx_defs.h'
  ],
  input: ['_generate_pyx.py', 'functions.json', '_add_newdocs.py'],
  command: [_generate_pyx, '-o', '@OUTDIR@'],
  install: true,
  install_dir: py3.get_install_dir() / 'scipy/special',
  install_tag: 'devel',
)

cython_special_pxd = fs.copyfile('cython_special.pxd')

# pyx -> c, pyx -> cpp generators, depending on copied pxi, pxd files.
uf_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_special_pxd])

uf_cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_special_pxd])

py3.extension_module('_ufuncs',
  [
    ufuncs_sources,
    uf_cython_gen.process(cython_special[0]),  # _ufuncs.pyx
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized] + scipy_import_dll_args,
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [
    lapack_dep,
    npymath_lib,
    np_dep,
    sf_error_state_dep,
  ],
  link_args: version_link_args,
  link_with: cdflib_lib,
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN',
)

ellint_files = [
  'ellint_carlson_wrap.hh',
  'ellint_carlson_cpp_lite/_rc.hh',
  'ellint_carlson_cpp_lite/_rd.hh',
  'ellint_carlson_cpp_lite/_rf.hh',
  'ellint_carlson_cpp_lite/_rg.hh',
  'ellint_carlson_cpp_lite/_rj.hh',
  'ellint_carlson_cpp_lite/ellint_argcheck.hh',
  'ellint_carlson_cpp_lite/ellint_arith_aux.hh',
  'ellint_carlson_cpp_lite/ellint_arithmetic.hh',
  'ellint_carlson_cpp_lite/ellint_carlson.hh',
  'ellint_carlson_cpp_lite/ellint_common.hh',
  'ellint_carlson_cpp_lite/ellint_typing.hh',
]

ellint_dep = declare_dependency(sources: ellint_files)

ufuncs_cxx_cpp_args = [
  cython_cpp_args,
  '-DBOOST_MATH_STANDALONE=1',
  '-DCYTHON_EXTERN_C=extern "C"',
]

py3.extension_module('_ufuncs_cxx',
  [ufuncs_cxx_sources,
    uf_cython_gen_cpp.process(cython_special[2]),  # _ufuncs_cxx.pyx
    ],
  cpp_args: ufuncs_cxx_cpp_args,
  include_directories: ['../_lib/boost_math/include', '../_lib',
                        '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [np_dep, ellint_dep, sf_error_state_dep],
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN',
)

py3.extension_module('_ellip_harm_2',
  [uf_cython_gen.process('_ellip_harm_2.pyx'), 'sf_error.cc'],
  c_args: cython_c_args,
  include_directories: ['../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [lapack_dep, np_dep, sf_error_state_dep],
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN',
)

py3.extension_module('cython_special',
  [
    uf_cython_gen.process('cython_special.pyx'),
    '_cosine.c',
    'special_wrappers.cpp',
    'sf_error.cc',
    'dd_real_wrappers.cpp'
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  include_directories: ['../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [np_dep, npymath_lib, sf_error_state_dep],
  link_with: cdflib_lib,
  install: true,
  subdir: 'scipy/special',
  install_rpath: '$ORIGIN',
)

py3.extension_module('_comb',
  cython_gen.process('_comb.pyx'),
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('_test_internal',
  [cython_gen.process('_test_internal.pyx'), 'dd_real_wrappers.cpp'],
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special',
  install_tag: 'tests'
)

py3.install_sources(
  '_test_internal.pyi',
  subdir: 'scipy/special',
  install_tag: 'tests'
)

# Must use `custom_target`, because `py3.install_sources` does not work with
# generated sources - see https://github.com/mesonbuild/meson/issues/7372
npz_files = [
  [
    '_data_boost',
    'tests/data/boost/assoc_legendre_p_ipp/assoc_legendre_p.txt',
    'boost',
    'boost.npz'
  ],
  [
    '_data_gsl',
    'tests/data/gsl/mathieu_ab.txt',
    'gsl',
    'gsl.npz'
  ],
  [
    '_data_local',
    'tests/data/local/ellipkm1.txt',
    'local',
    'local.npz'
  ],
]

foreach npz_file: npz_files
  custom_target(npz_file[0],
    output: npz_file[3],
    input: npz_file[1],
    command: [
      py3, '@CURRENT_SOURCE_DIR@/utils/makenpz.py',
      '--use-timestamp', npz_file[2], '-o', '@OUTDIR@'
    ],
    install: true,
    install_dir: py3.get_install_dir() / 'scipy/special/tests/data',
    install_tag: 'tests',
  )
endforeach


# Headers for special functions in `special` library (temporary name) are
# included in both build and install dirs for development purposes. Not public!
special_sources = [
  'special/config.h',
  'special/error.h',
  'special/evalpoly.h',
  'special/lambertw.h',
  'special/binom.h',
  'special/trig.h',
  'special/zlog1.h',
  'special/loggamma.h',
  'special/digamma.h',
  'special/wright_bessel.h',
  'special/hyp2f1.h',
  'special/tools.h',
  'special/iv_ratio.h'
]

special_cephes_sources = [
  'special/cephes/const.h',
  'special/cephes/polevl.h',
  'special/cephes/beta.h',
  'special/cephes/gamma.h',
  'special/cephes/trig.h',
  'special/cephes/zeta.h',
  'special/cephes/psi.h',
  'special/cephes/lanczos.h',
  'special/cephes/rgamma.h',
  'special/cephes/chbevl.h',
  'special/cephes/poch.h',
  'special/cephes/hyp2f1.h',
  'special/cephes/besselpoly.h',
  'special/cephes/i0.h',
  'special/cephes/i1.h',
  'special/cephes/scipy_iv.h',
  'special/cephes/j0.h',
  'special/cephes/j1.h',
  'special/cephes/jv.h',
  'special/cephes/cbrt.h',
  'special/cephes/airy.h',
  'special/cephes/k0.h',
  'special/cephes/k1.h',
  'special/cephes/kn.h',
  'special/cephes/ndtr.h',
  'special/cephes/igam.h',
  'special/cephes/unity.h',
  'special/cephes/chdtr.h',
  'special/cephes/igami.h',
  'special/cephes/hyperg.h',
  'special/cephes/expn.h'
]

py3.install_sources(special_sources, subdir: 'scipy/special/special')
py3.install_sources(special_cephes_sources, subdir: 'scipy/special/special/cephes')

python_sources = [
  '__init__.py',
  '_add_newdocs.py',
  '_basic.py',
  '_ellip_harm.py',
  '_lambertw.py',
  '_logsumexp.py',
  '_mptestutils.py',
  '_orthogonal.py',
  '_orthogonal.pyi',
  '_sf_error.py',
  '_spfun_stats.py',
  '_spherical_bessel.py',
  '_support_alternative_backends.py',
  '_testutils.py',
  '_ufuncs.pyi',
  'add_newdocs.py',
  'basic.py',
  'cython_special.pxd',
  'cython_special.pyi',
  'orthogonal.py',
  'sf_error.py',
  'specfun.py',
  'spfun_stats.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/special'
)

subdir('tests')
subdir('_precompute')
