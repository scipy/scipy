# TODO: 64-bit BLAS and LAPACK

_ufuncs_pxi_pxd_sources = custom_target('_ufuncs_pxi_pxd',
  output: [
    '__init__.py',
    '_agm.pxd',
    '_boxcox.pxd',
    '_cephes.pxd',
    '_complexstuff.pxd',
    '_convex_analysis.pxd',
    '_cunity.pxd',
    '_digamma.pxd',
    '_ellip_harm.pxd',
    '_ellip_harm_2.pxd',
    '_ellipk.pxd',
    '_evalpoly.pxd',
    '_exprel.pxd',
    '_factorial.pxd',
    '_hyp0f1.pxd',
    '_hyp2f1.pxd',
    '_hypergeometric.pxd',
    '_lambertw.pxd',
    '_legacy.pxd',
    '_loggamma.pxd',
    '_ndtri_exp.pxd',
    '_sici.pxd',
    '_spence.pxd',
    '_spherical_bessel.pxd',
    '_trig.pxd',
    '_wright_bessel.pxd',
    '_xlogy.pxd',
    'orthogonal_eval.pxd',
    'sf_error.pxd',
    'sph_harm.pxd',
    '_cython_special.pxi',
    '_cython_special_custom.pxi',
    '_ufuncs_extra_code.pxi',
    '_ufuncs_extra_code_common.pxi',
  ],
  input: [
    '__init__.py',
    '_agm.pxd',
    '_boxcox.pxd',
    '_cephes.pxd',
    '_complexstuff.pxd',
    '_convex_analysis.pxd',
    '_cunity.pxd',
    '_digamma.pxd',
    '_ellip_harm.pxd',
    '_ellip_harm_2.pxd',
    '_ellipk.pxd',
    '_evalpoly.pxd',
    '_exprel.pxd',
    '_factorial.pxd',
    '_hyp0f1.pxd',
    '_hyp2f1.pxd',
    '_hypergeometric.pxd',
    '_lambertw.pxd',
    '_legacy.pxd',
    '_loggamma.pxd',
    '_ndtri_exp.pxd',
    '_sici.pxd',
    '_spence.pxd',
    '_spherical_bessel.pxd',
    '_trig.pxd',
    '_wright_bessel.pxd',
    '_xlogy.pxd',
    'orthogonal_eval.pxd',
    'sf_error.pxd',
    'sph_harm.pxd',
    '_cython_special.pxi',
    '_cython_special_custom.pxi',
    '_ufuncs_extra_code.pxi',
    '_ufuncs_extra_code_common.pxi'
  ],
  command: [copier, '@INPUT@', '@OUTDIR@']
)

amos_sources = [
  'amos/dgamln.f',
  'amos/zabs.f',
  'amos/zacai.f',
  'amos/zacon.f',
  'amos/zairy.f',
  'amos/zasyi.f',
  'amos/zbesh.f',
  'amos/zbesi.f',
  'amos/zbesj.f',
  'amos/zbesk.f',
  'amos/zbesy.f',
  'amos/zbinu.f',
  'amos/zbiry.f',
  'amos/zbknu.f',
  'amos/zbuni.f',
  'amos/zbunk.f',
  'amos/zdiv.f',
  'amos/zexp.f',
  'amos/zkscl.f',
  'amos/zlog.f',
  'amos/zmlri.f',
  'amos/zmlt.f',
  'amos/zrati.f',
  'amos/zs1s2.f',
  'amos/zseri.f',
  'amos/zshch.f',
  'amos/zsqrt.f',
  'amos/zuchk.f',
  'amos/zunhj.f',
  'amos/zuni1.f',
  'amos/zuni2.f',
  'amos/zunik.f',
  'amos/zunk1.f',
  'amos/zunk2.f',
  'amos/zuoik.f',
  'amos/zwrsk.f'
]

cdflib_sources = [
  'cdflib/algdiv.f',
  'cdflib/alngam.f',
  'cdflib/alnrel.f',
  'cdflib/apser.f',
  'cdflib/basym.f',
  'cdflib/bcorr.f',
  'cdflib/betaln.f',
  'cdflib/bfrac.f',
  'cdflib/bgrat.f',
  'cdflib/bpser.f',
  'cdflib/bratio.f',
  'cdflib/brcmp1.f',
  'cdflib/brcomp.f',
  'cdflib/bup.f',
  'cdflib/cdfbet.f',
  'cdflib/cdfbin.f',
  'cdflib/cdfchi.f',
  'cdflib/cdfchn.f',
  'cdflib/cdff.f',
  'cdflib/cdffnc.f',
  'cdflib/cdfgam.f',
  'cdflib/cdfnbn.f',
  'cdflib/cdfnor.f',
  'cdflib/cdfpoi.f',
  'cdflib/cdft.f',
  'cdflib/cdftnc.f',
  'cdflib/cumbet.f',
  'cdflib/cumbin.f',
  'cdflib/cumchi.f',
  'cdflib/cumchn.f',
  'cdflib/cumf.f',
  'cdflib/cumfnc.f',
  'cdflib/cumgam.f',
  'cdflib/cumnbn.f',
  'cdflib/cumnor.f',
  'cdflib/cumpoi.f',
  'cdflib/cumt.f',
  'cdflib/cumtnc.f',
  'cdflib/devlpl.f',
  'cdflib/dinvnr.f',
  'cdflib/dinvr.f',
  'cdflib/dt1.f',
  'cdflib/dzror.f',
  'cdflib/erf.f',
  'cdflib/erfc1.f',
  'cdflib/esum.f',
  'cdflib/exparg.f',
  'cdflib/fpser.f',
  'cdflib/gam1.f',
  'cdflib/gaminv.f',
  'cdflib/gamln.f',
  'cdflib/gamln1.f',
  'cdflib/gamma_fort.f',
  'cdflib/grat1.f',
  'cdflib/gratio.f',
  'cdflib/gsumln.f',
  'cdflib/ipmpar.f',
  'cdflib/psi_fort.f',
  'cdflib/rcomp.f',
  'cdflib/rexp.f',
  'cdflib/rlog.f',
  'cdflib/rlog1.f',
  'cdflib/spmpar.f',
  'cdflib/stvaln.f'
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

ufuncs_sources = [
  '_cosine.c',
  'amos_wrappers.c',
  'cdf_wrappers.c',
  'specfun_wrappers.c',
  'sf_error.c'
]

ufuncs_cxx_sources = [
  '_faddeeva.cxx',
  '_wright.cxx',
  'ellint_carlson_wrap.cxx',
  'Faddeeva.cc',
  'sf_error.c',
  'wright.cc'
]


if is_windows
  use_math_defines = ['-D_USE_MATH_DEFINES']
else
  use_math_defines = []
endif

cephes_lib = static_library('cephes',
  cephes_sources,
  c_args: use_math_defines,
  include_directories: [inc_np, '../_lib', '../_build_utils/src'],
  dependencies: py3_dep
)

amos_lib = static_library('amos',
  amos_sources,
  fortran_args: '-Wno-maybe-uninitialized'
)

cdflib_lib = static_library('cdflib',
  cdflib_sources,
  fortran_args: [
    '-Wno-maybe-uninitialized',
    '-Wno-unused-label',
    '-Wno-intrinsic-shadow'
  ]
)

mach_lib = static_library('mach',
  mach_sources,
  fortran_args: ['-Wno-unused-dummy-argument', '-Wno-maybe-uninitialized']
)

specfun_lib = static_library('specfun',
  'specfun/specfun.f',
  fortran_args: '-Wno-maybe-uninitialized'
)

specfun_module = custom_target('specfun_module',
  output: '_specfunmodule.c',
  input: 'specfun.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_specfun',
  [specfun_module, fortranobject_c],
  c_args: numpy_nodepr_api,
  include_directories: [inc_np, inc_f2py],
  dependencies: py3_dep,
  link_with: specfun_lib,
  link_language: 'fortran',
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
  install_dir: py3.get_install_dir() + '/scipy/special'
)

# pyx -> c, pyx -> cpp generators, depending on copied pxi, pxd files.
uf_cython_gen = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@'],
  output : '@BASENAME@.c',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources])

uf_cython_gen_cpp = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@', '--cplus'],
  output : '@BASENAME@.cpp',
  depends : [_cython_tree, _ufuncs_pxi_pxd_sources])

py3.extension_module('_ufuncs',
  [
    ufuncs_sources,
    uf_cython_gen.process(cython_special[0]),  # _ufuncs.pyx
  ],
  c_args: [cython_c_args, use_math_defines, c_undefined_ok],
  include_directories: [inc_np, '../_lib', '../_build_utils/src'],
  dependencies: [
    py3_dep,
    lapack,
    npymath_lib,
    cython_tree
  ],
  link_with: [
    amos_lib,
    cdflib_lib,
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

py3.extension_module('_ufuncs_cxx',
  [ufuncs_cxx_sources,
    uf_cython_gen_cpp.process(cython_special[2]),  # _ufuncs_cxx.pyx
    ],
  cpp_args: [cython_c_args, use_math_defines],
  include_directories: [inc_np, '../_lib', '../_build_utils/src'],
  dependencies: [py3_dep, npymath_lib, ellint_dep],
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('_ellip_harm_2',
  [uf_cython_gen.process('_ellip_harm_2.pyx'), 'sf_error.c'],
  c_args: [cython_c_args, c_undefined_ok],
  include_directories: [inc_np, '../_lib', '../_build_utils/src'],
  dependencies: [py3_dep, lapack],
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('cython_special',
  [
    uf_cython_gen.process(cython_special[6]),  # cython_special.pyx
    '_cosine.c',
    'amos_wrappers.c',
    'cdf_wrappers.c',
    'specfun_wrappers.c',
    'sf_error.c'
  ],
  c_args: [cython_c_args, use_math_defines, c_undefined_ok],
  include_directories: [inc_np, '../_lib', '../_build_utils/src', 'cephes'],
  dependencies: [py3_dep, npymath_lib, lapack],
  link_with: [
    amos_lib,
    cdflib_lib,
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
  dependencies: py3_dep,
  install: true,
  subdir: 'scipy/special'
)

py3.extension_module('_test_round',
  cython_gen.process('_test_round.pyx'),
  include_directories: [inc_np, 'cephes', '../_lib', '../_build_utils/src'],
  dependencies: [py3_dep, npymath_lib],
  install: true,
  subdir: 'scipy/special'
)

# Must use `custom_target`, because `py3.install_sources` does not work with
# generated sources - see https://github.com/mesonbuild/meson/issues/7372
npz_files = [
  [
    '_data_boost',
    'tests/data/boost/acosh_data_ipp/acosh_data.txt',
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
    install_dir: py3.get_install_dir() + '/scipy/special/tests/data'
  )
endforeach


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
  '_test_round.pyi',
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
  pure: false,
  subdir: 'scipy/special'
)

subdir('tests')
subdir('_precompute')
