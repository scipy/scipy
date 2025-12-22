py3.extension_module('_vq',
  linalg_cython_gen.process('_vq.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/cluster/vq'
)

py3.install_sources([
    '__init__.py',
    '_vq_impl.py'
  ],
  subdir: 'scipy/cluster/vq'
)

subdir('tests')
