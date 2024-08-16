gen_sparsetools = find_program('../_generate_sparsetools.py')
_sparsetools_headers = custom_target('_sparsetools_headers',
  output: [
    'bsr_impl.h',
    'csc_impl.h',
    'csr_impl.h',
    'other_impl.h',
    'sparsetools_impl.h',
  ],
  input: '../_generate_sparsetools.py',
  command: [gen_sparsetools, '--no-force', '-o', '@OUTDIR@']
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
  dependencies: np_dep,
  include_directories: '../../_build_utils/src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse'
)
