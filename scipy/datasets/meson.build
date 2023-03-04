python_sources = [
  '__init__.py',
  '_fetchers.py',
  '_registry.py',
  '_download_all.py',
  '_utils.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/datasets'
)

if not get_option('disable-network-tests')
  subdir('tests')
endif
