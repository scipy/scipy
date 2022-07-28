python_sources = [
  '__init__.py',
]

_boost_gen_cpp = custom_target('_boost_gen_cpp',
  output: [
    'boost_ufuncs.cpp',
  ],
  input: '_gen.py',
  command: [py3, '@INPUT@', '-o', '@OUTDIR@'],
)

inc_boost = include_directories('..' / 'boost')

cpp_args = [
  #'-DBOOST_MATH_DOMAIN_ERROR_POLICY=ignore_error',
  #'-DBOOST_MATH_EVALUATION_ERROR_POLICY=user_error',
  #'-DBOOST_MATH_OVERFLOW_ERROR_POLICY=user_error'
]
if meson.get_compiler('cpp', native: true).sizeof('void*') > 4
  cpp_args += [numpy_nodepr_api]
endif

boostgen = py3.extension_module('_boostgen',
  _boost_gen_cpp,
  include_directories: [inc_np, inc_boost],
  cpp_args: [cpp_args],
  dependencies: [py3_dep, npymath_lib],
  install: true,
  subdir: 'scipy/_lib/boostgen'
)

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/_lib/boostgen'
)
