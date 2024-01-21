# TODO: 64-bit BLAS and LAPACK

_ufuncs_pxi_pxd_sources = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_agm.pxd'),
  fs.copyfile('_boxcox.pxd'),
  fs.copyfile('_cdflib.pxd'),
  fs.copyfile('_cdflib_wrappers.pxd'),
  fs.copyfile('_cephes.pxd'),
  fs.copyfile('_complexstuff.pxd'),
  fs.copyfile('_convex_analysis.pxd'),
  fs.copyfile('_cunity.pxd'),
  fs.copyfile('_digamma.pxd'),
  fs.copyfile('_ellip_harm.pxd'),
  fs.copyfile('_ellip_harm_2.pxd'),
  fs.copyfile('_ellipk.pxd'),
  fs.copyfile('_evalpoly.pxd'),
  fs.copyfile('_exprel.pxd'),
  fs.copyfile('_factorial.pxd'),
  fs.copyfile('_hyp0f1.pxd'),
  fs.copyfile('_hyp2f1.pxd'),
  fs.copyfile('_hypergeometric.pxd'),
  fs.copyfile('_legacy.pxd'),
  fs.copyfile('_loggamma.pxd'),
  fs.copyfile('_ndtri_exp.pxd'),
  fs.copyfile('_sici.pxd'),
  fs.copyfile('_spence.pxd'),
  fs.copyfile('_spherical_bessel.pxd'),
  fs.copyfile('_trig.pxd'),
  fs.copyfile('_wright_bessel.pxd'),
  fs.copyfile('_xlogy.pxd'),
  fs.copyfile('orthogonal_eval.pxd'),
  fs.copyfile('sf_error.pxd'),
  fs.copyfile('sph_harm.pxd'),
  fs.copyfile('_cython_special.pxi'),
  fs.copyfile('_cython_special_custom.pxi'),
  fs.copyfile('_ufuncs_extra_code.pxi'),
  fs.copyfile('_ufuncs_extra_code_common.pxi'),
]

cephes_sources = [
  'cephes/airy.c',
  'cephes/bdtr.c',
  'cephes/besselpoly.c',
  'cephes/beta.c',
  'cephes/btdtr.c',
  'cephes/cbrt.c',
  'cephes/chbevl.c',
  'cephes/chdtr.c',
  'cephes/const.c',
  'cephes/dawsn.c',
  'cephes/dd_real.c',
  'cephes/ellie.c',
  'cephes/ellik.c',
  'cephes/ellpe.c',
  'cephes/ellpj.c',
  'cephes/ellpk.c',
  'cephes/erfinv.c',
  'cephes/exp10.c',
  'cephes/exp2.c',
  'cephes/expn.c',
  'cephes/fdtr.c',
  'cephes/fresnl.c',
  'cephes/gamma.c',
  'cephes/gammasgn.c',
  'cephes/gdtr.c',
  'cephes/hyp2f1.c',
  'cephes/hyperg.c',
  'cephes/i0.c',
  'cephes/i1.c',
  'cephes/igam.c',
  'cephes/igami.c',
  'cephes/incbet.c',
  'cephes/incbi.c',
  'cephes/j0.c',
  'cephes/j1.c',
  'cephes/jv.c',
  'cephes/k0.c',
  'cephes/k1.c',
  'cephes/kn.c',
  'cephes/kolmogorov.c',
  'cephes/lanczos.c',
  'cephes/nbdtr.c',
  'cephes/ndtr.c',
  'cephes/ndtri.c',
  'cephes/owens_t.c',
  'cephes/pdtr.c',
  'cephes/poch.c',
  'cephes/psi.c',
  'cephes/rgamma.c',
  'cephes/round.c',
  'cephes/scipy_iv.c',
  'cephes/shichi.c',
  'cephes/sici.c',
  'cephes/sindg.c',
  'cephes/sinpi.c',
  'cephes/spence.c',
  'cephes/stdtr.c',
  'cephes/struve.c',
  'cephes/tandg.c',
  'cephes/tukey.c',
  'cephes/unity.c',
  'cephes/yn.c',
  'cephes/yv.c',
  'cephes/zeta.c',
  'cephes/zetac.c'
]

mach_sources = [
  'mach/d1mach.f',
  'mach/i1mach.f',
  'mach/xerror.f'
]

amos_lib = static_library('_amos',
  ['_amos.c'],
  c_args: use_math_defines,
  include_directories: ['../_lib', '../_build_utils/src'],
)

ufuncs_sources = [
  '_cosine.c',
  'scaled_exp1.c',
  'amos_wrappers.c',
  'specfun_wrappers.c',
  'sf_error.c'
]

ufuncs_cxx_sources = [
  '_faddeeva.cxx',
  '_wright.cxx',
  'ellint_carlson_wrap.cxx',
  'Faddeeva.cc',
  'sf_error.cc',
  'wright.cc'
]

cephes_lib = static_library('cephes',
  cephes_sources,
  c_args: use_math_defines,
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [py3_dep, np_dep],
)

mach_lib = static_library('mach',
  mach_sources,
  fortran_args: [_fflag_Wno_unused_dummy_argument, _fflag_Wno_maybe_uninitialized]
)

specfun_lib = static_library('specfun',
  'specfun/specfun.f',
  fortran_args: _fflag_Wno_maybe_uninitialized
)

specfun_module = custom_target('specfun_module',
  output: '_specfunmodule.c',
  input: 'specfun.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_specfun',
  specfun_module,
  link_args: version_link_args,
  dependencies: [fortranobject_dep],
  link_with: specfun_lib,
  link_language: 'fortran',
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('_cdflib',
  cython_gen.process('_cdflib.pyx'),
  c_args: cython_c_args,
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
cython_special = custom_target('cython_special',
  output: [
    '_ufuncs.pyx',
    '_ufuncs_defs.h',
    '_ufuncs_cxx.pyx',
    '_ufuncs_cxx.pxd',
    '_ufuncs_cxx_defs.h',
    '_ufuncs.pyi',
    'cython_special.pyx',
    'cython_special.pxd'
  ],
  input: ['_generate_pyx.py', 'functions.json', '_add_newdocs.py'],
  command: [py3, '@INPUT0@', '-o', '@OUTDIR@'],
  install: true,
  install_dir: py3.get_install_dir() / 'scipy/special',
  install_tag: 'devel',
)

# Only needed to establish a dependency; see comments in linalg/meson.build
# for cython_blas_pxd
cython_special_pxd = custom_target('_dummy_cython_special.pxd',
  output: '_dummy_cython_special.pxd',
  input: cython_special[7],
  command: [py3, '-c', '"@INPUT@"'],
)

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
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  include_directories: ['../_lib', '../_build_utils/src'],
  dependencies: [
    lapack,
    npymath_lib,
    np_dep,
  ],
  link_args: version_link_args,
  link_with: [
    amos_lib,
    cephes_lib,
    mach_lib,
    specfun_lib
  ],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/special'
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
  '-DSP_SPECFUN_ERROR', # For error handling in special functions.
]

py3.extension_module('_ufuncs_cxx',
  [ufuncs_cxx_sources,
    uf_cython_gen_cpp.process(cython_special[2]),  # _ufuncs_cxx.pyx
    ],
  cpp_args: ufuncs_cxx_cpp_args,
  include_directories: ['../_lib/boost_math/include', '../_lib',
                        '../_build_utils/src'],
  link_args: version_link_args,
    link_with: [cephes_lib],
  dependencies: [np_dep, ellint_dep],
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('_ellip_harm_2',
  [uf_cython_gen.process('_ellip_harm_2.pyx'), 'sf_error.c'],
  c_args: cython_c_args,
  include_directories: ['../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [lapack, np_dep],
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('cython_special',
  [
    uf_cython_gen.process(cython_special[6]),  # cython_special.pyx
    '_cosine.c',
    'scaled_exp1.c',
    'amos_wrappers.c',
    'specfun_wrappers.c',
    'sf_error.c'
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  include_directories: ['../_lib', '../_build_utils/src', 'cephes'],
  link_args: version_link_args,
  dependencies: [np_dep, npymath_lib, lapack],
  link_with: [
    amos_lib,
    cephes_lib,
    mach_lib,
    specfun_lib
  ],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/special'
)

py3.extension_module('_comb',
  cython_gen.process('_comb.pyx'),
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special'
)

dd_dep = declare_dependency(sources: ['cephes/dd_real.c'])

py3.extension_module('_test_internal',
  cython_gen.process('_test_internal.pyx'),
  include_directories: ['cephes', '../_lib', '../_build_utils/src'],
  dependencies: [np_dep, dd_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/special'
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
    install_tag: 'python-runtime',
  )
endforeach


# Headers for special functions in `special` library (temporary name) are
# included in both build and install dirs for development purposes. Not public!
special_sources = [
  'special/config.h',
  'special/error.h',
  'special/evalpoly.h',
  'special/lambertw.h',
  'special/binom.h'
]

special_cephes_sources = [
  'special/cephes/const.h',
  'special/cephes/polevl.h',
  'special/cephes/beta.h',
  'special/cephes/gamma.h'
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
  '_test_internal.pyi',
  '_testutils.py',
  'add_newdocs.py',
  'basic.py',
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
