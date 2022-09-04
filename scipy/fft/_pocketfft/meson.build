thread_dep = dependency('threads', required: false)

# mingw-w64 does not implement pthread_pfork, needed by pocket_fft
win_gcc = is_windows and meson.get_compiler('cpp').get_id() == 'gcc'

pocketfft_threads = []
fft_deps = []
if not win_gcc
  if thread_dep.found()
    fft_deps += [thread_dep]
    if not is_windows
      # `pthreads` probably not installed, and native threading is always
      # available. Not easy to distinguish this better, Meson builtin
      # functionality for that in progress (see comment on gh-16957).
      # The code on `pocketfft_hdronly.h` will include `<threads>` anyway.
      pocketfft_threads += ['-DPOCKETFFT_PTHREADS']
    endif
  endif
else
  # Freezes using threading for mingw-w64 gcc:
  # https://github.com/mreineck/pocketfft/issues/1
  pocketfft_threads += ['-DPOCKETFFT_NO_MULTITHREADING']
endif

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
  'LICENSE.md',
  'realtransforms.py'
]


py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy/fft/_pocketfft'
)

subdir('tests')
