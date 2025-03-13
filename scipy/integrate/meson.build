mach_src = [
  'mach/d1mach.f',
  'mach/xerror.f'
]

lsoda_src = [
  'odepack/blkdta000.f',
  'odepack/bnorm.f',
  'odepack/cfode.f',
  'odepack/ewset.f',
  'odepack/fnorm.f',
  'odepack/intdy.f',
  'odepack/lsoda.f',
  'odepack/prja.f',
  'odepack/solsy.f',
  'odepack/srcma.f',
  'odepack/stoda.f',
  'odepack/vmnorm.f',
  'odepack/xerrwv.f',
  'odepack/xsetf.f',
  'odepack/xsetun.f'
]

vode_src = [
  'odepack/vode.f',
  'odepack/zvode.f'
]

dop_src = [
  'dop/dop853.f',
  'dop/dopri5.f'
]

quadpack_test_src = [
  'tests/_test_multivariate.c'
]


if use_ilp64
  _fflag_intsize = _fflag_ilp64
else
  _fflag_intsize = _fflag_lp64
endif


mach_lib = static_library('mach_lib',
  mach_src,
  gnu_symbol_visibility: 'hidden',
  fortran_args: [fortran_ignore_warnings, _fflag_intsize],
)


lsoda_lib = static_library('lsoda_lib',
  lsoda_src,
  fortran_args: [fortran_ignore_warnings, _fflag_intsize],
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

vode_lib = static_library('vode_lib',
  vode_src,
  fortran_args: [fortran_ignore_warnings, _fflag_intsize],
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

dop_lib = static_library('dop_lib',
  dop_src,
  fortran_args: [fortran_ignore_warnings, _fflag_intsize],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('_quadpack',
  ['__quadpack.h', '__quadpack.c'],
  dependencies: [np_dep, ccallback_dep, lapack_ilp64],
  install: true,
  subdir: 'scipy/integrate'
)


_c_args_lsoda = []
if use_ilp64
  _c_args_lsoda += ['-DHAVE_BLAS_ILP64']
endif

py3.extension_module('_odepack',
  '_odepackmodule.c',
  fortran_args: [fortran_ignore_warnings, _fflag_intsize],
  c_args: _c_args_lsoda,
  link_with: [lsoda_lib, mach_lib],
  link_args: version_link_args,
  dependencies: [lapack_ilp64, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)


if use_ilp64
   # generator only accepts strings, not files
   f2c_map_file = f2py_ilp64_opts[1]
   extra_arg = f2py_ilp64_opts[0] + '=' + fs.parent(f2c_map_file) / fs.name(f2c_map_file)
  _vodemodule_obj = f2py_gen.process('vode.pyf', extra_args: extra_arg)
  _lsodamodule_obj = f2py_gen.process('lsoda.pyf', extra_args: extra_arg)
  _dopmodule_obj = f2py_gen.process('dop.pyf', extra_args: extra_arg)
  _testodeintmodule_obj = f2py_gen.process('tests/test_odeint_banded.pyf', extra_args: extra_arg)
else
  _vodemodule_obj = f2py_gen.process('vode.pyf')
  _lsodamodule_obj = f2py_gen.process('lsoda.pyf')
  _dopmodule_obj = f2py_gen.process('dop.pyf')
  _testodeintmodule_obj = f2py_gen.process('tests/test_odeint_banded.pyf')
endif


py3.extension_module('_vode',
  _vodemodule_obj,
  link_with: [vode_lib],
  fortran_args: [_fflag_intsize],
  c_args: [Wno_unused_variable] + _c_args_lsoda,
  link_args: version_link_args,
  dependencies: [lapack_ilp64, blas_ilp64, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_lsoda',
  _lsodamodule_obj,
  link_with: [lsoda_lib, mach_lib],
  fortran_args: [_fflag_intsize],
  c_args: [Wno_unused_variable] + _c_args_lsoda,
  dependencies: [lapack_ilp64, fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_dop',
  _dopmodule_obj,
  link_with: [dop_lib],
  fortran_args: [_fflag_intsize],
  c_args: [Wno_unused_variable] + _c_args_lsoda,
  dependencies: [fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_test_multivariate',
  [quadpack_test_src],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate',
  install_tag: 'tests'
)

py3.extension_module('_test_odeint_banded',
  ['tests/banded5x5.f', _testodeintmodule_obj],
  link_with: [lsoda_lib, mach_lib],
  fortran_args: [_fflag_intsize, _fflag_Wno_unused_dummy_argument],
  c_args: _c_args_lsoda,
  link_args: version_link_args,
  dependencies: [lapack_ilp64, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate',
  install_tag: 'tests'
)

subdir('_ivp')
subdir('_rules')
subdir('tests')

py3.install_sources([
    '__init__.py',
    '_bvp.py',
    '_lebedev.py',
    '_ode.py',
    '_quad_vec.py',
    '_quadrature.py',
    '_odepack_py.py',
    '_quadpack_py.py',
    '_tanhsinh.py',
    '_cubature.py',
    'dop.py',
    'lsoda.py',
    'odepack.py',
    'quadpack.py',
    'vode.py',
  ],
  subdir: 'scipy/integrate'
)
