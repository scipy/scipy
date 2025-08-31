propack_sources = [
    'PROPACK/src/common.c',
    'PROPACK/src/gemm_overwrite.c',
    'PROPACK/src/getu0.c',
    'PROPACK/src/gs.c',
    'PROPACK/src/lanbpro.c',
    'PROPACK/src/lansvd_irl.c',
    'PROPACK/src/lansvd.c',
    'PROPACK/src/ritzvec.c'
]

propack_incdir_public = include_directories('PROPACK/include')
propack_incdir_private = include_directories('PROPACK/src/include')

propack_lib = static_library('propack_lib',
  propack_sources,
  include_directories: [propack_incdir_public, propack_incdir_private],
  dependencies: [m_dep, lapack_lp64_dep],
)

py3.extension_module('_propack',
  '_propackmodule.c',
  link_with: propack_lib,
  dependencies: [np_dep, ccallback_dep],
  link_args: [version_link_args],
  install: true,
  subdir: 'scipy/sparse/linalg/',
)
