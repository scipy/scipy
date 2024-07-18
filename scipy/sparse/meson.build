_csparsetools_pyx = custom_target('_csparsetools_pyx',
  output: '_csparsetools.pyx',
  input: '_csparsetools.pyx.in',
  command: [tempita, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_csparsetools',
  cython_gen.process(_csparsetools_pyx),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse'
)

python_sources = [
  '__init__.py',
  '_base.py',
  '_bsr.py',
  '_compressed.py',
  '_construct.py',
  '_coo.py',
  '_csc.py',
  '_csr.py',
  '_data.py',
  '_dia.py',
  '_dok.py',
  '_extract.py',
  '_index.py',
  '_lil.py',
  '_matrix_io.py',
  '_matrix.py',
  '_spfuncs.py',
  '_sputils.py',
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
  subdir: 'scipy/sparse'
)

subdir('sparsetools')
subdir('csgraph')
subdir('linalg')
subdir('tests')
