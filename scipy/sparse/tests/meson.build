python_sources = [
  '__init__.py',
  'test_arithmetic1d.py',
  'test_array_api.py',
  'test_base.py',
  'test_coo.py',
  'test_common1d.py',
  'test_construct.py',
  'test_csc.py',
  'test_csr.py',
  'test_dok.py',
  'test_extract.py',
  'test_matrix_io.py',
  'test_minmax1d.py',
  'test_sparsetools.py',
  'test_spfuncs.py',
  'test_sputils.py'
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/sparse/tests',
  install_tag: 'tests'
)

data_sources = [
  'data/csc_py2.npz',
  'data/csc_py3.npz'
]

py3.install_sources(
  data_sources,
  subdir: 'scipy/sparse/tests/data',
  install_tag: 'tests'
)
