project('biteopt',
  version : '2024.6',
  meson_version: '>= 1.5.0',
)

fs = import('fs')
if not fs.exists('biteopt/README.md')
  error('Missing the `biteopt` submodule! Run `git submodule update --init` to fix this.')
endif

biteopt_dep = declare_dependency(
  include_directories: include_directories('biteopt')
)
