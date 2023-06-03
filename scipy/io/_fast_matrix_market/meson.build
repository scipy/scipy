subdir('fast_matrix_market/dependencies/fast_float')
subdir('fast_matrix_market/dependencies/ryu')

py3.extension_module('_core',
  ['src/fast_matrix_market/_core.cpp'],
  override_options: ['cpp_std=c++17'],
  cpp_args: [numpy_nodepr_api, '-DFMM_FROM_CHARS_INT_SUPPORTED', '-DFMM_TO_CHARS_INT_SUPPORTED', '-DFMM_USE_FAST_FLOAT', '-DFMM_USE_RYU'],
  include_directories: ['fast_matrix_market/include'],
  dependencies: [np_dep, pybind11_dep, fast_float_dep, ryu_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/io/_fast_matrix_market'
)

py3.install_sources([
    'src/fast_matrix_market/__init__.py'
  ],
  subdir: 'scipy/io/_fast_matrix_market'
)
