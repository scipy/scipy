cpp_args = [
  '-DBOOST_MATH_DOMAIN_ERROR_POLICY=ignore_error',
  '-DBOOST_MATH_EVALUATION_ERROR_POLICY=user_error',
  '-DBOOST_MATH_OVERFLOW_ERROR_POLICY=user_error'
]

if meson.get_compiler('cpp', native: true).sizeof('void*') > 4
  cpp_args += ['-DBOOST_MATH_PROMOTE_DOUBLE_POLICY=false',
                numpy_nodepr_api]
endif

include = include_directories('include')
inc_boost = include_directories('..' / '..' / '_lib' / 'boost')

hypergeom_ufunc = py3.extension_module('hypergeom_ufunc',
  hypergeom_ufunc_pyx,
  include_directories: [include, inc_np, inc_boost],
  cpp_args: [cpp_args, cython_cpp_args],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/stats/_boost'
)

nbinom_ufunc = py3.extension_module('nbinom_ufunc',
  nbinom_ufunc_pyx,
  include_directories: [include, inc_np, inc_boost],
  cpp_args: [cpp_args, cython_cpp_args],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/stats/_boost'
)

beta_ufunc = py3.extension_module('beta_ufunc',
  beta_ufunc_pyx,
  include_directories: [include, inc_np, inc_boost],
  cpp_args: [cpp_args, cython_cpp_args],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/stats/_boost'
)

binom_ufunc = py3.extension_module('binom_ufunc',
  binom_ufunc_pyx,
  include_directories: [include, inc_np, inc_boost],
  cpp_args: [cpp_args, cython_cpp_args],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/stats/_boost'
)

ncf_ufunc = py3.extension_module('ncf_ufunc',
  ncf_ufunc_pyx,
  include_directories: [include, inc_np, inc_boost],
  cpp_args: [cpp_args, cython_cpp_args],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/stats/_boost'
)

py3.install_sources([
    '__init__.py'
  ],
  pure: false,
  subdir: 'scipy/stats/_boost'
)
