py3.extension_module('_sparsetools',
  [
    'bsr.cxx',
    'csc.cxx',
    'csr.cxx',
    'other.cxx',
    'sparsetools.cxx',
  ],
  dependencies: [np_dep, pybind11_dep],
  include_directories: '../../_build_utils/src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse'
)
