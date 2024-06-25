python_sources = [
  '__init__.py',
  'dots.png',
  'test_c_api.py',
  'test_datatypes.py',
  'test_filters.py',
  'test_fourier.py',
  'test_interpolation.py',
  'test_measurements.py',
  'test_morphology.py',
  'test_ni_support.py',
  'test_splines.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/ndimage/tests',
  install_tag: 'tests'
)

py3.install_sources([
    'data/label_inputs.txt',
    'data/label_results.txt',
    'data/label_strels.txt',
  ],
  subdir: 'scipy/ndimage/tests/data',
  install_tag: 'tests'
)
