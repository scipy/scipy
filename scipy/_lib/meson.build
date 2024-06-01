fs = import('fs')
if not fs.exists('boost_math/README.md')
  error('Missing the `boost` submodule! Run `git submodule update --init` to fix this.')
endif
if not fs.exists('highs/README.md')
  error('Missing the `highs` submodule! Run `git submodule update --init` to fix this.')
endif
if not fs.exists('unuran/README.md')
  error('Missing the `unuran` submodule! Run `git submodule update --init` to fix this.')
endif
if not fs.exists('array_api_compat/README.md')
  error('Missing the `array_api_compat` submodule! Run `git submodule update --init` to fix this.')
endif
if not fs.exists('pocketfft/README.md')
  error('Missing the `pocketfft` submodule! Run `git submodule update --init` to fix this.')
endif
if not fs.exists('cobyqa/README.rst')
  error('Missing the `cobyqa` submodule! Run `git submodule update --init` to fix this.')
endif

_lib_pxd = [
  fs.copyfile('__init__.py'),
  fs.copyfile('_ccallback_c.pxd'),
  fs.copyfile('ccallback.pxd'),
  fs.copyfile('messagestream.pxd'),
]

# Cython pyx -> c generator with _lib_pxd dependency
lib_cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : [_cython_tree, _lib_pxd])

py3.extension_module('_ccallback_c',
  lib_cython_gen.process('_ccallback_c.pyx'),
  c_args: [cython_c_args, Wno_discarded_qualifiers],
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib'
)

py3.extension_module('_test_ccallback',
  'src/_test_ccallback.c',
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib',
  install_tag: 'tests'
)

py3.extension_module('_fpumode',
  '_fpumode.c',
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib'
)

py3.extension_module('_test_deprecation_call',
  cython_gen.process('_test_deprecation_call.pyx'),
  c_args: cython_c_args,
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib',
  install_tag: 'tests'
)

py3.extension_module('_test_deprecation_def',
  cython_gen.process('_test_deprecation_def.pyx'),
  c_args: cython_c_args,
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib',
  install_tag: 'tests'
)

# May be easier as a compile flag, but use a config header file to stay
# in sync with what is done in setup.py
# TODO: the `prefix` here is recommended in
#       https://mesonbuild.com/Compiler-properties.html#does-a-function-exist
#       and seems needed, but why does it fail to detect it for the conda
#       compiler? Does it not exist, or is it hidden elsewhere?

conf_memstream = configuration_data()

if meson.get_compiler('c').has_function('open_memstream',
                                        prefix: '#include <stdio.h>')
  conf_memstream.set('has_openmemstream', '1')
else
  conf_memstream.set('has_openmemstream', '0')
endif

configure_file(
  input: 'src/messagestream_config.h.in',
  output: 'messagestream_config.h',
  configuration: conf_memstream
)

py3.extension_module('messagestream',
  lib_cython_gen.process('messagestream.pyx'),
  c_args: cython_c_args,
  include_directories: 'src',
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib'
)

python_sources = [
  '__init__.py',
  '_array_api.py',
  '_bunch.py',
  '_ccallback.py',
  '_disjoint_set.py',
  '_docscrape.py',
  '_elementwise_iterative_method.py',
  '_finite_differences.py',
  '_gcutils.py',
  '_pep440.py',
  '_testutils.py',
  '_threadsafety.py',
  '_tmpdirs.py',
  '_util.py',
  'decorator.py',
  'deprecation.py',
  'doccer.py',
  'uarray.py',
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/_lib'
)

# `array_api_compat` install to simplify import path;
# should be updated whenever new files are added to `array_api_compat`

py3.install_sources(
  [
    'array_api_compat/array_api_compat/__init__.py',
    'array_api_compat/array_api_compat/_internal.py',
  ],
  subdir: 'scipy/_lib/array_api_compat',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/common/__init__.py',
    'array_api_compat/array_api_compat/common/_aliases.py',
    'array_api_compat/array_api_compat/common/_fft.py',
    'array_api_compat/array_api_compat/common/_helpers.py',
    'array_api_compat/array_api_compat/common/_linalg.py',
    'array_api_compat/array_api_compat/common/_typing.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/common',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/cupy/__init__.py',
    'array_api_compat/array_api_compat/cupy/_aliases.py',
    'array_api_compat/array_api_compat/cupy/_typing.py',
    'array_api_compat/array_api_compat/cupy/fft.py',
    'array_api_compat/array_api_compat/cupy/linalg.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/cupy',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/dask/__init__.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/dask/',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/dask/array/__init__.py',
    'array_api_compat/array_api_compat/dask/array/_aliases.py',
    'array_api_compat/array_api_compat/dask/array/linalg.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/dask/array',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/numpy/__init__.py',
    'array_api_compat/array_api_compat/numpy/_aliases.py',
    'array_api_compat/array_api_compat/numpy/_typing.py',
    'array_api_compat/array_api_compat/numpy/fft.py',
    'array_api_compat/array_api_compat/numpy/linalg.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/numpy',
)

py3.install_sources(
  [
    'array_api_compat/array_api_compat/torch/__init__.py',
    'array_api_compat/array_api_compat/torch/_aliases.py',
    'array_api_compat/array_api_compat/torch/fft.py',
    'array_api_compat/array_api_compat/torch/linalg.py',
  ],
  subdir: 'scipy/_lib/array_api_compat/torch',
)

py3.install_sources(
  [
    'cobyqa/cobyqa/__init__.py',
    'cobyqa/cobyqa/framework.py',
    'cobyqa/cobyqa/main.py',
    'cobyqa/cobyqa/models.py',
    'cobyqa/cobyqa/problem.py',
    'cobyqa/cobyqa/settings.py',
  ],
  subdir: 'scipy/_lib/cobyqa',
)

py3.install_sources(
  [
    'cobyqa/cobyqa/subsolvers/__init__.py',
    'cobyqa/cobyqa/subsolvers/geometry.py',
    'cobyqa/cobyqa/subsolvers/optim.py',
  ],
  subdir: 'scipy/_lib/cobyqa/subsolvers',
)

py3.install_sources(
  [
    'cobyqa/cobyqa/utils/__init__.py',
    'cobyqa/cobyqa/utils/exceptions.py',
    'cobyqa/cobyqa/utils/math.py',
    'cobyqa/cobyqa/utils/versions.py',
  ],
  subdir: 'scipy/_lib/cobyqa/utils',
)

subdir('_uarray')
subdir('tests')
