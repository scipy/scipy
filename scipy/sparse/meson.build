_csparsetools_pyx = custom_target('_csparsetools_pyx',
  output: '_csparsetools.pyx',
  input: '_csparsetools.pyx.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_csparsetools',
  _csparsetools_pyx,
  include_directories: inc_np,
  dependencies: py3_dep,
  install: true,
  subdir: 'scipy/sparse'
)

python_sources = [
  '__init__.py',
  '_index.py',
  '_matrix_io.py',
  'base.py',
  'bsr.py',
  'compressed.py',
  'construct.py',
  'coo.py',
  'csc.py',
  'csr.py',
  'data.py',
  'dia.py',
  'dok.py',
  'extract.py',
  'lil.py',
  'sparsetools.py',
  'spfuncs.py',
  'sputils.py'
]


py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/sparse'
)

subdir('sparsetools')
subdir('csgraph')
subdir('linalg')
subdir('tests')
