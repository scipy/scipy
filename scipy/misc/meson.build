python_sources = [
  '__init__.py',
  '_common.py',
  'common.py',
  'doccer.py'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/misc'
)

py3.install_sources([
    'ascent.dat',
    'ecg.dat',
    'face.dat'
  ],
  pure: false,
  subdir: 'scipy/misc'
)

subdir('tests')
