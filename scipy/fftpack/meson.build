py3.extension_module('convolve',
  cython_gen.process('convolve.pyx'),
  c_args: cython_c_args,
  include_directories: inc_np,
  dependencies: py3_dep,
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
  pure: false,             # Will be installed next to binaries
  subdir: 'scipy/fftpack'  # Folder relative to site-packages to install to
)

subdir('tests')
