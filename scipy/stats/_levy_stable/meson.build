py3.install_sources([
    '__init__.py',
  ],
  subdir: 'scipy/stats/_levy_stable'
)

rustmod = import('rust')

cargo_ws = rustmod.workspace()
pkg = cargo_ws.package()

py3.extension_module('levyst',
  'src/levyst.rs',
  rust_dependency_map : pkg.rust_dependency_map(),
  rust_args: pkg.rust_args(),
  dependencies : pkg.dependencies(),
  install: true,
  subdir: 'scipy/stats/_levy_stable',
)
