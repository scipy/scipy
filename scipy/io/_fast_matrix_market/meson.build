subdir('fast_matrix_market/dependencies/fast_float')
subdir('fast_matrix_market/dependencies/ryu')

if meson.get_compiler('cpp').get_id() == 'msvc'
  fmm_extra_args = ['/bigobj']
else
  fmm_extra_args = []
endif

py3.extension_module('_core',
  [
    '_core.hpp',
    '_core.cpp',
    '_core_read_array.cpp',
    '_core_read_coo.cpp',
    '_core_write_array.cpp',
    '_core_write_coo.cpp',
    '_core_write_csc.cpp',
  ],
  override_options: ['cpp_std=c++17'],
  cpp_args: [numpy_nodepr_api, '-DFMM_FROM_CHARS_INT_SUPPORTED', '-DFMM_TO_CHARS_INT_SUPPORTED', '-DFMM_USE_FAST_FLOAT', '-DFMM_USE_RYU'] + fmm_extra_args,
  include_directories: ['fast_matrix_market/include'],
  dependencies: [np_dep, pybind11_dep, fast_float_dep, ryu_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/io/_fast_matrix_market'
)

py3.install_sources([
    '__init__.py'
  ],
  subdir: 'scipy/io/_fast_matrix_market'
)
