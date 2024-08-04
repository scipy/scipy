python_sources = [
  '__init__.py',
  '_common.py',
  'common.py',
  'doccer.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/misc'
)

py3.install_sources([
    'ascent.dat',
    'ecg.dat',
    'face.dat'
  ],
  subdir: 'scipy/misc'
)

subdir('tests')
