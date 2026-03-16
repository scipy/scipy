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

py3.extension_module('_propack',
  [
    propack_sources,
    '_propackmodule.c'
  ],
  include_directories: ['PROPACK/include', 'PROPACK/src/include'],
  dependencies: [np_dep, ccallback_dep, lapack_lp64_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse/linalg/',
)
