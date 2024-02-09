python_sources = [
  '__init__.py',
  '_codata.py',
  '_constants.py',
  'codata.py',
  'constants.py'
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/constants'
)

py3.install_sources([
  codata_constants_2002.txt,
  codata_constants_2006.txt,
  codata_constants_2010.txt,
  codata_constants_2014.txt,
  codata_constants_2018.txt
  ],
subdir: 'scipy/constants'
)


subdir('tests')
