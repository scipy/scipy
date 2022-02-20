_sparsetools_headers = custom_target('_sparsetools_headers',
  output: [
    'bsr_impl.h',
    'csc_impl.h',
    'csr_impl.h',
    'other_impl.h',
    'sparsetools_impl.h',
  ],
  input: '../_generate_sparsetools.py',
  command: [py3, '@INPUT@', '--no-force', '-o', '@OUTDIR@']
)

py3.extension_module('_sparsetools',
  [
    'bsr.cxx',
    'csc.cxx',
    'csr.cxx',
    'other.cxx',
    'sparsetools.cxx',
    _sparsetools_headers,
  ],
  cpp_args: numpy_nodepr_api,
  include_directories: inc_np,
  dependencies: py3_dep,
  install: true,
  subdir: 'scipy/sparse'
)
