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

mach_lib = static_library('mach_lib',
  mach_src,
  fortran_args: fortran_ignore_warnings,
  gnu_symbol_visibility: 'hidden',
)

lsoda_lib = static_library('lsoda_lib',
  lsoda_src,
  fortran_args: fortran_ignore_warnings,
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

vode_lib = static_library('vode_lib',
  vode_src,
  fortran_args: fortran_ignore_warnings,
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

dop_lib = static_library('dop_lib',
  dop_src,
  fortran_args: fortran_ignore_warnings,
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('_quadpack',
  ['__quadpack.h', '__quadpack.c'],
  include_directories: ['../_lib/src'],
  dependencies: np_dep,
  install: true,
  subdir: 'scipy/integrate'
)

py3.extension_module('_odepack',
  '_odepackmodule.c',
  fortran_args: fortran_ignore_warnings,
  link_with: [lsoda_lib, mach_lib],
  link_args: version_link_args,
  dependencies: [lapack_dep, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_vode',
  f2py_gen.process('vode.pyf'),
  link_with: [vode_lib],
  c_args: [Wno_unused_variable],
  link_args: version_link_args,
  dependencies: [lapack_dep, blas_dep, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_lsoda',
  f2py_gen.process('lsoda.pyf'),
  link_with: [lsoda_lib, mach_lib],
  c_args: [Wno_unused_variable],
  dependencies: [lapack_dep, fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_dop',
  f2py_gen.process('dop.pyf'),
  link_with: [dop_lib],
  c_args: [Wno_unused_variable],
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
  ['tests/banded5x5.f', f2py_gen.process('tests/test_odeint_banded.pyf')],
  link_with: [lsoda_lib, mach_lib],
  fortran_args: _fflag_Wno_unused_dummy_argument,
  link_args: version_link_args,
  dependencies: [lapack_dep, fortranobject_dep],
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
