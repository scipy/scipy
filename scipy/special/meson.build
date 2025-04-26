# TODO: 64-bit BLAS and LAPACK

_ufuncs_pxi_pxd_sources = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_agm.pxd'),
  fs.copyfile('_boxcox.pxd'),
  fs.copyfile('_cdflib_wrappers.pxd'),
  fs.copyfile('_complexstuff.pxd'),
  fs.copyfile('_convex_analysis.pxd'),
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
  fs.copyfile('orthogonal_eval.pxd'),
  fs.copyfile('sf_error.pxd'),
  fs.copyfile('_ufuncs_extra_code.pxi'),
  fs.copyfile('_ufuncs_extra_code_common.pxi'),
]

ufuncs_sources = [
  '_cosine.c',
  'xsf_wrappers.cpp',
  'sf_error.cc',
  'dd_real_wrappers.cpp'
]

ufuncs_cxx_sources = [
  '_wright.cxx',
  'ellint_carlson_wrap.cxx',
  'sf_error.cc',
  'wright.cc'
]

cdflib_lib = static_library('cdflib',
  'cdflib.c',
  include_directories: ['../_lib', '../_build_utils/src'],
  gnu_symbol_visibility: 'hidden',
)

# The -Wno-psabi arg is used to silence a warning that the ABI of
# passing structure with complex float member changed in GCC 4.4,
# which is not relevant for scipy.special.
ufuncs_cpp_args = ['-DSP_SPECFUN_ERROR']
if meson.get_compiler('cpp').get_id() == 'gcc'
  ufuncs_cpp_args += ['-Wno-psabi']
endif

py3.extension_module('_special_ufuncs',
  ['_special_ufuncs.cpp', '_special_ufuncs_docs.cpp', 'sf_error.cc'],
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  dependencies: [xsf_dep, np_dep],
  link_args: version_link_args,
  cpp_args: ufuncs_cpp_args,
  install: true,
  subdir: 'scipy/special',
)

py3.extension_module('_gufuncs',
  ['_gufuncs.cpp', '_gufuncs_docs.cpp', 'sf_error.cc'],
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  dependencies: [xsf_dep, np_dep],
  link_args: version_link_args,
  cpp_args: ufuncs_cpp_args,
  install: true,
  subdir: 'scipy/special',
)

py3.extension_module('_specfun',
  [cython_gen_cpp.process('_specfun.pyx')],
  dependencies: [xsf_dep, np_dep],
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
_cython_tree += fs.copyfile('__init__.pxd')

# pyx -> c, pyx -> cpp generators, depending on copied pxi, pxd files.
uf_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_special_pxd, cython_special[0]])

uf_cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources, cython_special_pxd, cython_special[0]])

py3.extension_module('_ufuncs',
  [
    ufuncs_sources,
    uf_cython_gen.process(cython_special[0]),  # _ufuncs.pyx
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  cpp_args: ['-DSP_SPECFUN_ERROR'],
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  dependencies: [
    xsf_dep,
    lapack_dep,
    npymath_lib,
    np_dep,
  ],
  link_args: version_link_args,
  link_with: cdflib_lib,
  install: true,
  subdir: 'scipy/special',
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
  '-DSP_SPECFUN_ERROR'
]

py3.extension_module('_ufuncs_cxx',
  [ufuncs_cxx_sources,
    uf_cython_gen_cpp.process(cython_special[2]),  # _ufuncs_cxx.pyx
    ],
  cpp_args: ufuncs_cxx_cpp_args,
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [boost_math_dep, xsf_dep, np_dep, ellint_dep],
  install: true,
  subdir: 'scipy/special',
)

py3.extension_module('_ellip_harm_2',
  [uf_cython_gen.process('_ellip_harm_2.pyx'), 'sf_error.cc'],
  c_args: [cython_c_args],
  cpp_args: ['-DSP_SPECFUN_ERROR'],
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [xsf_dep, lapack_dep, np_dep],
  install: true,
  subdir: 'scipy/special',
)

py3.extension_module('cython_special',
  [
    uf_cython_gen.process('cython_special.pyx'),
    '_cosine.c',
    'xsf_wrappers.cpp',
    'sf_error.cc',
    'dd_real_wrappers.cpp'
  ],
  c_args: [cython_c_args, Wno_maybe_uninitialized],
  cpp_args: ['-DSP_SPECFUN_ERROR'],
  include_directories: ['..', '../_lib', '../_build_utils/src'],
  link_args: version_link_args,
  dependencies: [xsf_dep, np_dep, npymath_lib, lapack_dep],
  link_with: cdflib_lib,
  install: true,
  subdir: 'scipy/special',
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
  dependencies: [xsf_dep, np_dep],
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


python_sources = [
  '__init__.py',
  '__init__.pxd',
  '_add_newdocs.py',
  '_basic.py',
  '_ellip_harm.py',
  '_input_validation.py',
  '_lambertw.py',
  '_logsumexp.py',
  '_mptestutils.py',
  '_multiufuncs.py',
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
