py3.extension_module('_sigtools',
  [
    '_firfilter.c',
    '_sigtoolsmodule.c',
    '_medianfilter.c',
    tempita_gen.process('_correlate_nd.c.in'),
    tempita_gen.process('_lfilter.c.in'),
  ],
  dependencies: np_dep,
  link_args: version_link_args,
  include_directories: ['../_build_utils/src'],
  install: true,
  subdir: 'scipy/signal'
)

if use_pythran
  py3.extension_module('_max_len_seq_inner',
    pythran_gen.process('_max_len_seq_inner.py'),
    cpp_args: [cpp_args_pythran, _cpp_Wno_unused_local_typedefs],
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/signal'
  )

  py3.extension_module('_spectral',
    pythran_gen.process('_spectral.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/signal'
  )
else
  py3.extension_module('_max_len_seq_inner',
    [cython_gen.process('_max_len_seq_inner.pyx')],
    c_args: cython_c_args,
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/signal'
  )

  py3.extension_module('_spectral',
    [cython_gen.process('_spectral.pyx')],
    c_args: cython_c_args,
    dependencies: np_dep,
    link_args: version_link_args,
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
    dependencies: np_dep,
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/signal'
  )
endforeach

py3.extension_module('_spline',
  [
    '_splinemodule.c',
    tempita_gen.process('_bspline_util.c.in'),
  ],
  dependencies: np_dep,
  link_args: version_link_args,
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
    '_peak_finding.py',
    '_savitzky_golay.py',
    '_short_time_fft.py',
    '_signaltools.py',
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
  subdir: 'scipy/signal'
)

subdir('windows')
subdir('tests')
