py3.extension_module('_interpnd',
  spt_cython_gen.process('_interpnd.pyx'),
  c_args: [Wno_discarded_qualifiers, cython_c_args],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

py3.extension_module('_ppoly',
  linalg_cython_gen.process('_ppoly.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

py3.extension_module('_rgi_cython',
  linalg_cython_gen.process('_rgi_cython.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  install: true,
  subdir: 'scipy/interpolate'
)

__fitpack_lib = static_library('__fitpack',
    ['src/__fitpack.h', 'src/__fitpack.cc'],
    dependencies:[lapack_dep, np_dep, py3_dep],
)

__fitpack_dep = declare_dependency(
    link_with: __fitpack_lib,
)


py3.extension_module('_dierckx',
    ['src/_dierckxmodule.cc'],
    include_directories: 'src/',
    dependencies: [np_dep, __fitpack_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/interpolate'
)

py3.extension_module('_fitpack',
  ['src/_fitpackmodule.c', 'src/dfitpack.c'],
  include_directories: 'src/',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

if use_pythran
  # cannot use pythran_gen: generator does not support renaming the source -> extension
  # there must be a better way
  _rbf_src = custom_target('_rbf_src',
     output: '_rbfinterp_pythran.cpp',
     input: '_rbfinterp_pythran_src.py',
     command: [pythran, '-E', '@INPUT@', '-o', '@OUTPUT@']
   )

  py3.extension_module('_rbfinterp_pythran',
    _rbf_src,
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/interpolate'
  )
else
  # XXX need to rename _rbfinterp_pythran_src.py in the source tree into _rbfinterp_pythran.py
  py3.install_sources(
    ['_rbfinterp_pythran.py'],
    subdir: 'scipy/interpolate'
  )
endif

py3.install_sources([
    '__init__.py',
    '_bsplines.py',
    '_bary_rational.py',
    '_cubic.py',
    '_fitpack_impl.py',
    '_fitpack_py.py',
    '_fitpack2.py',
    '_fitpack_repro.py',
    '_interpolate.py',
    '_ndgriddata.py',
    '_pade.py',
    '_polyint.py',
    '_rbf.py',
    '_rbfinterp.py',
    '_rbfinterp_np.py',
    '_rbfinterp_xp.py',
    '_rbfinterp_pythran_src.py',
    '_rgi.py',
    '_ndbspline.py',
    'dfitpack.py',
    'fitpack.py',
    'fitpack2.py',
    'interpolate.py',
    'ndgriddata.py',
    'interpnd.py',
    'polyint.py',
    'rbf.py'
  ],
  subdir: 'scipy/interpolate'
)

subdir('tests')
