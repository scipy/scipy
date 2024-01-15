pocketfft_threads = []
fft_deps = []
if is_windows
  if is_mingw
    # mingw-w64 does not implement pthread_pfork, needed by pocket_fft
    # Disable threading completely, because of freezes using threading for
    # mingw-w64 gcc: https://github.com/mreineck/pocketfft/issues/1
    pocketfft_threads += ['-DPOCKETFFT_NO_MULTITHREADING']
  else
    if thread_dep.found()
      # Use native Windows threading for MSVC/Clang. `pthreads` is probably not
      # installed, and native threading is always available. It is not easy to
      # distinguish this better, Meson builtin functionality for that is in
      # progress (see comment on gh-16957).  The code in `pocketfft_hdronly.h`
      # will include `<threads>` anyway.
      fft_deps += [thread_dep]
      pocketfft_threads += []
    endif
  endif
else
  if thread_dep.found()
    fft_deps += [thread_dep]
    pocketfft_threads += ['-DPOCKETFFT_PTHREADS']
  endif
endif
pocketfft_args = pocketfft_threads + ['-DPOCKETFFT_CACHE_SIZE=16']

incdir = include_directories('../../_lib/pocketfft')

py3.extension_module('pypocketfft',
  'pypocketfft.cxx',
  cpp_args: pocketfft_args,
  dependencies: [fft_deps, pybind11_dep],
  link_args: version_link_args,
  include_directories: incdir,
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
  subdir: 'scipy/fft/_pocketfft'
)

subdir('tests')
