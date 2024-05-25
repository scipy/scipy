py3.extension_module('_nd_image',
  [
    'src/nd_image.c',
    'src/ni_filters.c',
    'src/ni_fourier.c',
    'src/ni_interpolation.c',
    'src/ni_measure.c',
    'src/ni_morphology.c',
    'src/ni_splines.c',
    'src/ni_support.c'
  ],
  include_directories: ['../_lib/src', '../_build_utils/src'],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/ndimage'
)

py3.extension_module('_ni_label',
  cython_gen.process('src/_ni_label.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/ndimage'
)

py3.extension_module('_ctest',
  'src/_ctest.c',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/ndimage',
  install_tag: 'tests'
)

py3.extension_module('_cytest',
  cython_gen.process('src/_cytest.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/ndimage',
  install_tag: 'tests'
)


python_sources = [
  '__init__.py',
  '_filters.py',
  '_fourier.py',
  '_interpolation.py',
  '_measurements.py',
  '_morphology.py',
  '_ni_docstrings.py',
  '_ni_support.py',
  'filters.py',
  'fourier.py',
  'interpolation.py',
  'measurements.py',
  'morphology.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/ndimage'
)

subdir('tests')
