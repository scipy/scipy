python_sources = [
  '__init__.py',
  'test_basic.py',
  'test_helper.py',
  'test_import.py',
  'test_pseudo_diffs.py',
  'test_real_transforms.py'
]

test_sources = [
  'fftw_double_ref.npz',
  'fftw_longdouble_ref.npz',
  'fftw_single_ref.npz',
  'test.npz'
]

py3.install_sources(
  [python_sources, test_sources],
  pure: false,
  subdir: 'scipy/fftpack/tests'
)
