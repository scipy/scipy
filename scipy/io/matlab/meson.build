streams = py3.extension_module('streams',
  ['streams.pyx'],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/io/matlab'
)

mio_utils = py3.extension_module('mio_utils',
  ['mio_utils.pyx'],
#  c_args: numpy_nodepr_api, TODO: Uncomment once NumPy API usage is fixed
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/io/matlab'
)

mio5_utils = py3.extension_module('mio5_utils',
  ['mio5_utils.pyx'],
#  c_args: numpy_nodepr_api, TODO: Uncomment once NumPy API usage is fixed
  include_directories: [incdir_numpy],
  dependencies: [py3_dep],
  install: true,
  subdir: 'scipy/io/matlab'
)

py3.install_sources([
    '__init__.py',
    'byteordercodes.py',
    'mio.py',
    'mio4.py',
    'mio5.py',
    'mio5_params.py',
    'miobase.py'
  ],
  pure: false,
  subdir: 'scipy/io/matlab'
)

subdir('tests')
