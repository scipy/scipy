# Setup the highs library
highs_proj = subproject('highs',
                      default_options : ['default_library=static',
                                        'with_pybind11=True'])
highs_dep = highs_proj.get_variable('highs_dep')
highspy_py = highs_proj.get_variable('highspy_py')
highspy_cpp = highs_proj.get_variable('highspy_cpp')
highsoptions_cpp = highs_proj.get_variable('highsoptions_cpp')

pyb11_dep = [
  py3.dependency(),
  dependency('pybind11')
]

py3.extension_module(
  '_highs',
  sources : highspy_cpp,
  dependencies: [pyb11_dep, highs_dep],
  c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
  cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
  subdir: 'highspy',
  install: true,
)


py3.extension_module(
  '_highs_options',
  sources : highsoptions_cpp,
  dependencies: [pyb11_dep, highs_dep],
  c_args: [Wno_unused_variable, Wno_unused_but_set_variable],
  cpp_args: [_cpp_Wno_unused_variable, _cpp_Wno_unused_but_set_variable],
  subdir: 'highspy',
  install: true,
)

py3.install_sources(
  highspy_py,
  subdir: 'highspy',
)

py3.install_sources([
  '__init__.py',
  '_highs_wrapper.py',
],
  subdir: 'scipy/optimize/_highs'
)
