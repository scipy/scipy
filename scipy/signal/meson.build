correlate_nd_c = custom_target('_correlate_nd',
  output: '_correlate_nd.c',
  input: '_correlate_nd.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

lfilter_c = custom_target('_lfilter',
  output: '_lfilter.c',
  input: '_lfilter.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

sigtools = py3.extension_module('_sigtools',
  [
    '_sigtoolsmodule.c',
    '_firfilter.c',
    '_sigtools.h',
    '_medianfilter.c',
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
    cpp_args: [cpp_args_pythran, _cpp_Wno_unused_local_typedefs],
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
    cpp_args: cpp_args_pythran,
    include_directories: [incdir_pythran, incdir_numpy],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )
else
  _max_len_seq_inner = py3.extension_module('_max_len_seq_inner',
    [cython_gen.process('_max_len_seq_inner.pyx')],
    c_args: cython_c_args,
    include_directories: [inc_np],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )

  _spectral = py3.extension_module('_spectral',
    [cython_gen.process('_spectral.pyx')],
    c_args: cython_c_args,
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
    cython_gen.process(pyx_file[1]),
    c_args: cython_c_args,
    include_directories: [inc_np],
    dependencies: [py3_dep],
    install: true,
    subdir: 'scipy/signal'
  )
endforeach

bspline_util = custom_target('_bspline_util',
  output: '_bspline_util.c',
  input: '_bspline_util.c.in',
  command: [py3, tempita, '@INPUT@', '-o', '@OUTDIR@']
)

spline = py3.extension_module('_spline',
  ['_splinemodule.c', bspline_util],
  c_args: numpy_nodepr_api,
  include_directories: [inc_np],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/signal'
)

py3.install_sources([
    '__init__.py',
    '_arraytools.py',
    '_bsplines.py',
    '_czt.py',
    '_filter_design.py',
    '_fir_filter_design.py',
    '_lti_conversion.py',
    '_ltisys.py',
    '_max_len_seq.py',
    '_max_len_seq_inner.py',
    '_peak_finding.py',
    '_savitzky_golay.py',
    '_signaltools.py',
    '_spectral.py',
    '_spectral_py.py',
    '_upfirdn.py',
    '_waveforms.py',
    '_wavelets.py',
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
