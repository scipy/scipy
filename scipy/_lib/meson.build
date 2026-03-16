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

ccallback_dep = declare_dependency(include_directories: ['..', 'src'])

py3.extension_module('_ccallback_c',
  lib_cython_gen.process('_ccallback_c.pyx'),
  c_args: [cython_c_args, Wno_discarded_qualifiers],
  dependencies: ccallback_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/_lib'
)

py3.extension_module('_test_ccallback',
  'src/_test_ccallback.c',
  dependencies: ccallback_dep,
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
                                        prefix: '#define _POSIX_C_SOURCE 200809L\n#include <stdio.h>')
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
  '_array_api_docs_tables.py',
  '_array_api_no_0d.py',
  '_array_api_override.py',
  '_bunch.py',
  '_ccallback.py',
  '_disjoint_set.py',
  '_docscrape.py',
  '_elementwise_iterative_method.py',
  '_gcutils.py',
  '_public_api.py',
  '_sparse.py',
  '_testutils.py',
  '_util.py',
  'deprecation.py',
  'doccer.py',
  'uarray.py',
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/_lib'
)

subdir('_uarray')
subdir('tests')
