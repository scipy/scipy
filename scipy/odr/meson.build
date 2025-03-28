odrpack = static_library('odrpack',
  [
    'odrpack/d_lpk.f',
    'odrpack/d_mprec.f',
    'odrpack/d_odr.f',
    'odrpack/dlunoc.f'
  ],
  fortran_args: _fflag_Wno_conversion,  # silence "conversion from REAL(8) to INTEGER(4)"
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('__odrpack',
  '__odrpack.c',
  link_with: odrpack,
  link_args: version_link_args,
  dependencies: [blas_dep, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/odr'
)


python_sources = [
  '__init__.py',
  '_add_newdocs.py',
  '_models.py',
  '_odrpack.py',
  'models.py',
  'odrpack.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/odr'
)

subdir('tests')
