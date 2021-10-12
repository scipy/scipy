correlate_nd_c = custom_target('correlate_nd',
  output: 'correlate_nd.c',
  input: 'correlate_nd.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

lfilter_c = custom_target('lfilter',
  output: 'lfilter.c',
  input: 'lfilter.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

sigtools = py3.extension_module('sigtools',
  [
    'sigtoolsmodule.c',
    'firfilter.c',
    'sigtools.h',
    'medianfilter.c',
    lfilter_c,
    correlate_nd_c
  ],
  c_args: numpy_nodepr_api,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/signal'
)

if use_pythran
  _max_len_seq_inner = custom_target('_max_len_seq_inner',
    output: ['_max_len_seq_inner.cpp'],
    input: '_max_len_seq_inner.py',
    command: [pythran, '-E', '@INPUT@', '-o', '@OUTDIR@/_max_len_seq_inner.cpp']
  )

  _max_len_seq_inner = py3.extension_module('_max_len_seq_inner',
    [_max_len_seq_inner],
    cpp_args: [
      '-Wno-unused-function', '-Wno-unused-variable',
      '-Wno-deprecated-declarations', '-Wno-unused-local-typedefs',
      '-Wno-cpp', '-Wno-int-in-bool-context'
    ] + cpp_args_pythran,
    include_directories: [incdir_pythran, incdir_numpy],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )

  _spectral = custom_target('_spectral',
    output: ['_spectral.cpp'],
    input: '_spectral.py',
    command: [pythran, '-E', '@INPUT@', '-o', '@OUTDIR@/_spectral.cpp']
  )

  _spectral = py3.extension_module('_spectral',
    [_spectral],
    cpp_args: [
      '-Wno-unused-function', '-Wno-unused-variable',
      '-Wno-deprecated-declarations',
      '-Wno-cpp', '-Wno-int-in-bool-context'
    ] + cpp_args_pythran,
    include_directories: [incdir_pythran, incdir_numpy],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )
else
  _max_len_seq_inner = py3.extension_module('_max_len_seq_inner',
    ['_max_len_seq_inner.pyx'],
    include_directories: [inc_np],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )

  _spectral = py3.extension_module('_spectral',
    ['_spectral.pyx'],
    include_directories: [inc_np],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )
endif

pyx_files = [
  ['_peak_finding_utils', '_peak_finding_utils.pyx'],
  ['_sosfilt', '_sosfilt.pyx'],
  ['_upfirdn_apply', '_upfirdn_apply.pyx']
]

foreach pyx_file: pyx_files
  py3.extension_module(pyx_file[0],
    pyx_file[1],
    include_directories: [inc_np],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )
endforeach

bspline_util = custom_target('bspline_util',
  output: 'bspline_util.c',
  input: 'bspline_util.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

spline = py3.extension_module('_spline',
  ['splinemodule.c', bspline_util],
  c_args: numpy_nodepr_api,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/signal'
)

py3.install_sources([
    '__init__.py',
    '_arraytools.py',
    '_max_len_seq.py',
    '_max_len_seq_inner.py',
    '_peak_finding.py',
    '_savitzky_golay.py',
    '_spectral.py',
    '_upfirdn.py',
    'bsplines.py',
    'filter_design.py',
    'fir_filter_design.py',
    'lti_conversion.py',
    'ltisys.py',
    'signaltools.py',
    'spectral.py',
    'spline.py',
    'waveforms.py',
    'wavelets.py'
  ],
  pure: false,
  subdir: 'scipy/signal'
)

subdir('windows')
subdir('tests')
