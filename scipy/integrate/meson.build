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
  fortran_args: fortran_ignore_warnings,
  gnu_symbol_visibility: 'hidden',
)

quadpack_lib = static_library('quadpack_lib',
  quadpack_src,
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
  ['__quadpack.h', '_quadpackmodule.c'],
  fortran_args: fortran_ignore_warnings,
  include_directories: ['../_lib/src'],
  link_with: [quadpack_lib, mach_lib],
  link_args: version_link_args,
  dependencies: [lapack_dep, np_dep],
  install: true,
  link_language: 'fortran',
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
    'dop.py',
    'lsoda.py',
    'odepack.py',
    'quadpack.py',
    'vode.py',
  ],
  subdir: 'scipy/integrate'
)
