project('boost-math',
  version : '1.88.0',
  meson_version: '>= 1.5.0',
)

fs = import('fs')
if not fs.exists('math/README.md')
  error('Missing the `boost_math` submodule! Run `git submodule update --init` to fix this.')
endif

boost_math_dep = declare_dependency(include_directories: 'math/include')
