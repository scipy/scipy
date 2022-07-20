python_sources = [
  '__init__.py',
  '_backend.py',
  '_basic.py',
  '_debug_backends.py',
  '_fftlog.py',
  '_fftlog_multimethods.py',
  '_helper.py',
  '_realtransforms.py'
]


py3.install_sources(
  python_sources,
  pure: false,         # Will be installed next to binaries
  subdir: 'scipy/fft'  # Folder relative to site-packages to install to
)

subdir('_pocketfft')
subdir('tests')
