py3.extension_module('convolve',
  cython_gen.process('convolve.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/fftpack'
)


python_sources = [
  '__init__.py',
  '_basic.py',
  '_helper.py',
  '_pseudo_diffs.py',
  '_realtransforms.py',
  'basic.py',
  'helper.py',
  'pseudo_diffs.py',
  'realtransforms.py'
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/fftpack'
)

subdir('tests')
