subdir('fast_matrix_market/dependencies/fast_float')
subdir('fast_matrix_market/dependencies/ryu')

if meson.get_compiler('cpp').get_id() == 'msvc'
  fmm_extra_args = ['/bigobj']
else
  fmm_extra_args = []
endif

# macOS gained <charconv> support with macOS 10.15.
# As of writing SciPy targets macOS 10.9 on x86 (and 12 on arm64).
# Enable use of <charconv> only if not on macOS.
if host_machine.system() != 'darwin'
  fmm_extra_args += [
    '-DFMM_FROM_CHARS_INT_SUPPORTED',
    '-DFMM_TO_CHARS_INT_SUPPORTED'
    ]
endif

py3.extension_module('_fmm_core',
  [
    'src/_fmm_core.hpp',
    'src/_fmm_core.cpp',
    'src/_fmm_core_read_array.cpp',
    'src/_fmm_core_read_coo.cpp',
    'src/_fmm_core_write_array.cpp',
    'src/_fmm_core_write_coo_32.cpp',
    'src/_fmm_core_write_coo_64.cpp',
  ],
  override_options: ['cpp_std=c++17'],
  cpp_args: [
    numpy_nodepr_api,
    '-DFMM_SCIPY_PRUNE',
    '-DFMM_USE_FAST_FLOAT',
    '-DFMM_USE_RYU'
    ] + fmm_extra_args,
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
