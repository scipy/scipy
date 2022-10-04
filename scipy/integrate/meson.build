mach_src = [
  'mach/d1mach.f',
  'mach/xerror.f'
]

quadpack_src = [
  'quadpack/dqag.f',
  'quadpack/dqage.f',
  'quadpack/dqagi.f',
  'quadpack/dqagie.f',
  'quadpack/dqagp.f',
  'quadpack/dqagpe.f',
  'quadpack/dqags.f',
  'quadpack/dqagse.f',
  'quadpack/dqawc.f',
  'quadpack/dqawce.f',
  'quadpack/dqawf.f',
  'quadpack/dqawfe.f',
  'quadpack/dqawo.f',
  'quadpack/dqawoe.f',
  'quadpack/dqaws.f',
  'quadpack/dqawse.f',
  'quadpack/dqc25c.f',
  'quadpack/dqc25f.f',
  'quadpack/dqc25s.f',
  'quadpack/dqcheb.f',
  'quadpack/dqelg.f',
  'quadpack/dqk15.f',
  'quadpack/dqk15i.f',
  'quadpack/dqk15w.f',
  'quadpack/dqk21.f',
  'quadpack/dqk31.f',
  'quadpack/dqk41.f',
  'quadpack/dqk51.f',
  'quadpack/dqk61.f',
  'quadpack/dqmomo.f',
  'quadpack/dqng.f',
  'quadpack/dqpsrt.f',
  'quadpack/dqwgtc.f',
  'quadpack/dqwgtf.f',
  'quadpack/dqwgts.f'
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
  fortran_args: fortran_ignore_warnings
)

quadpack_lib = static_library('quadpack_lib',
  quadpack_src,
  fortran_args: fortran_ignore_warnings
)

lsoda_lib = static_library('lsoda_lib',
  lsoda_src,
  fortran_args: fortran_ignore_warnings
)

vode_lib = static_library('vode_lib',
  vode_src,
  fortran_args: fortran_ignore_warnings
)

dop_lib = static_library('dop_lib',
  dop_src,
  fortran_args: fortran_ignore_warnings
)

py3.extension_module('_quadpack',
  ['__quadpack.h', '_quadpackmodule.c'],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  include_directories: [inc_np, '../_lib/src'],
  link_with: [quadpack_lib, mach_lib],
  dependencies: [py3_dep, lapack],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_odepack',
  '_odepackmodule.c',
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  link_with: [lsoda_lib, mach_lib],
  dependencies: [py3_dep, lapack],
  include_directories: [inc_np],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

vode_module = custom_target('vode_module',
  output: ['_vode-f2pywrappers.f', '_vodemodule.c'],
  input: 'vode.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_vode',
  [vode_module, fortranobject_c],
  link_with: [vode_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  include_directories: [inc_np, inc_f2py],
  dependencies: [py3_dep, lapack],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

lsoda_module = custom_target('lsoda_module',
  output: ['_lsoda-f2pywrappers.f', '_lsodamodule.c'],
  input: 'lsoda.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_lsoda',
  [lsoda_module, fortranobject_c],
  link_with: [lsoda_lib, mach_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  dependencies: [py3_dep, lapack],
  include_directories: [inc_np, inc_f2py],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

_dop_module = custom_target('_dop_module',
  output: ['_dop-f2pywrappers.f', '_dopmodule.c'],
  input: 'dop.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_dop',
  [_dop_module, fortranobject_c],
  link_with: [dop_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  dependencies: [py3_dep],
  include_directories: [inc_np, inc_f2py],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_test_multivariate',
  [quadpack_test_src],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/integrate'
)

_test_odeint_banded_module = custom_target('_test_odeint_banded_module',
  output: ['_test_odeint_bandedmodule.c', '_test_odeint_banded-f2pywrappers.f'],
  input: 'tests/banded5x5.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_test_odeint_banded',
  ['tests/banded5x5.f', fortranobject_c, _test_odeint_banded_module],
  c_args: numpy_nodepr_api,
  link_with: [lsoda_lib, mach_lib],
  fortran_args: _fflag_Wno_unused_dummy_argument,
  include_directories: [inc_np, inc_f2py],
  dependencies: [py3_dep, lapack],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

subdir('_ivp')
subdir('tests')

py3.install_sources([
    '__init__.py',
    '_bvp.py',
    '_ode.py',
    '_quad_vec.py',
    '_quadrature.py',
    '_odepack_py.py',
    '_quadpack_py.py',
    'dop.py',
    'lsoda.py',
    'odepack.py',
    'quadpack.py',
    'setup.py',
    'vode.py',
  ],
  pure: false,
  subdir: 'scipy/integrate'
)
