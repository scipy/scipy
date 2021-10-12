thread_dep = dependency('threads', required: false)

if thread_dep.found()
  pocketfft_threads = '-DPOCKETFFT_PTHREADS'
endif

#TODO: add the equivalent of set_cxx_flags_hook
py3.extension_module('pypocketfft',
  'pypocketfft.cxx',
  cpp_args: pocketfft_threads,
  include_directories: inc_pybind11,
  dependencies: [py3_dep, thread_dep],
  gnu_symbol_visibility: 'hidden',
  install: true,
  subdir: 'scipy/fft/_pocketfft',
)


python_sources = [
  '__init__.py',
  'basic.py',
  'helper.py',
  'realtransforms.py'
]


py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/fft/_pocketfft'
)

subdir('tests')
