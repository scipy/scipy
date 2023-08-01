python_sources = [
  '__init__.py',
  '_backend.py',
  '_basic.py',
  '_basic_uarray.py',
  '_debug_backends.py',
  '_fftlog.py',
  '_fftlog_np.py',
  '_fftlog_multimethods.py',
  '_helper.py',
  '_realtransforms.py',
  '_realtransforms_uarray.py'
]


py3.install_sources(
  python_sources,
  subdir: 'scipy/fft'
)

subdir('_pocketfft')
subdir('tests')
