thread_dep = dependency('threads', required: false)

# mingw-w64 does not implement pthread_pfork, needed by pocket_fft
win_gcc = is_windows and meson.get_compiler('cpp').get_id() == 'gcc'

pocketfft_threads = []
fft_deps = []
if not win_gcc
  if thread_dep.found()
    pocketfft_threads += ['-DPOCKETFFT_PTHREADS']
    fft_deps += [thread_dep]
  endif
else
  # Freezes using threading for mingw-w64 gcc:
  # https://github.com/mreineck/pocketfft/issues/1
  pocketfft_threads += ['-DPOCKETFFT_NO_MULTITHREADING']
endif

#TODO: add the equivalent of set_cxx_flags_hook
py3.extension_module('pypocketfft',
  'pypocketfft.cxx',
  cpp_args: pocketfft_threads,
  include_directories: inc_pybind11,
  dependencies: fft_deps,
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
