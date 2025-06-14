propack_sources = [
    '_propack/src/common.c',
    '_propack/src/gemm_overwrite.c',
    '_propack/src/getu0.c',
    '_propack/src/gs.c',
    '_propack/src/lanbpro.c',
    '_propack/src/lansvd_irl.c',
    '_propack/src/lansvd.c',
    '_propack/src/ritzvec.c'
]

propack_incdir_pub = include_directories('_propack/include')
propack_incdir_priv = include_directories('_propack/src/include')

propack_lib = static_library('propack',
  propack_sources,
  include_directories: [propack_incdir_pub, propack_incdir_priv],
  dependencies: [m_dep, lapack_dep, blas_dep],
)

py3.extension_module('_propacklib',
  '_propackmodule.c',
  link_with: propack_lib,
  dependencies: [np_dep],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/sparse/linalg/_propack',
)

py3.install_sources([
    '__init__.py',
    '_propack.py',
  ],
  subdir: 'scipy/sparse/linalg/_eigen/_propack'
)

subdir('tests')
