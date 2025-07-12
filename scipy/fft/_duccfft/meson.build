duccfft_threads = []
fft_deps = []
if is_windows
  if is_mingw
    # mingw-w64 does not implement pthread_pfork, needed by ducc_fft
    # Disable threading completely, because of freezes using threading for
    # mingw-w64 gcc: https://github.com/mreineck/pocketfft/issues/1
    duccfft_threads += ['-DDUCC0_NO_LOWLEVEL_THREADING']
  else
    if thread_dep.found()
      # Use native Windows threading for MSVC/Clang. `pthreads` is probably not
      # installed, and native threading is always available. It is not easy to
      # distinguish this better, Meson builtin functionality for that is in
      # progress (see comment on gh-16957).  The ducc FFT sources
      # will include `<threads>` anyway.
      fft_deps += [thread_dep]
      duccfft_threads += []
    endif
  endif
else
  if thread_dep.found()
    fft_deps += [thread_dep]
  endif
endif
duccfft_args = duccfft_threads

incdir = include_directories('../../_lib/ducc/src')

py3.extension_module('pyduccfft',
  '../../_lib/ducc/src/ducc0/infra/string_utils.cc',
  '../../_lib/ducc/src/ducc0/infra/threading.cc',
  '../../_lib/ducc/src/ducc0/infra/mav.cc',
  'pyduccfft.cxx',
  cpp_args: duccfft_args,
  dependencies: [fft_deps, pybind11_dep],
  link_args: version_link_args,
  include_directories: incdir,
  install: true,
  subdir: 'scipy/fft/_duccfft',
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
  subdir: 'scipy/fft/_duccfft'
)

subdir('tests')
