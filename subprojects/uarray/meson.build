project(
  'uarray',
  'cpp',
  version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
  license: 'BSD-3-Clause',
  license_files: ['LICENSE'],
  default_options: 'cpp_std=c++17',
)

fs = import('fs')
py = import('python').find_installation(pure: false)

# Sources variable definitions #

# NOTE: downstream projects may rely on these variable names when vendoring,
# do not change them without consulting with downstream projects
python_sources = {
  'uarray': files(
    'src/uarray/__init__.py',
    'src/uarray/_backend.py',
    'src/uarray/_typing.pyi',
    'src/uarray/_typing.pyi',
    'src/uarray/_version.pyi',
  )
}

if fs.exists('src/uarray/_version.py')
  python_sources = python_sources + {
    'uarray': python_sources['uarray'] + files('src/uarray/_version.py')
  }
endif

cpp_sources = files('src/_uarray_dispatch.cxx')
include_dirs = include_directories('src')

# Installation

if not meson.is_subproject()
  py.extension_module('_uarray',
    cpp_sources,
    cpp_args: meson.get_compiler('cpp').get_supported_arguments('-Wno-terminate'),
    include_directories: include_dirs,
    install: true,
    subdir: 'uarray',
  )
endif

subdir('src/uarray')
