# Setup the highs library
fs = import('fs')
if not fs.exists('../../../subprojects/highs/README.md')
  error('Missing the `highs` submodule! Run `git submodule update --init` to fix this.')
endif
highs_proj = subproject('highs',
                      default_options : ['default_library=static',
                                        'use_zlib=disabled'])
highs_dep = highs_proj.get_variable('highs_dep')
highspy_cpp = highs_proj.get_variable('highspy_cpp')
highsoptions_cpp = highs_proj.get_variable('highsoptions_cpp')

scipy_highspy_dep = [
  py3_dep,
  pybind11_dep,
  highs_dep,
  thread_dep,
  atomic_dep,
]

py3.extension_module(
  '_highs',
  sources : highspy_cpp,
  dependencies: scipy_highspy_dep,
  c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
  cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
  link_args: version_link_args,
  subdir: 'scipy/optimize/_highs/highspy',
  install: true,
)


py3.extension_module(
  '_highs_options',
  sources : highsoptions_cpp,
  dependencies: scipy_highspy_dep,
  c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
  cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
  link_args: version_link_args,
  subdir: 'scipy/optimize/_highs/highspy',
  install: true,
)

py3.install_sources([
  '__init__.py',
  '_highs_wrapper.py',
],
  subdir: 'scipy/optimize/_highs'
)
