fitpack_src = [
  'fitpack/bispeu.f',
  'fitpack/bispev.f',
  'fitpack/clocur.f',
  'fitpack/cocosp.f',
  'fitpack/concon.f',
  'fitpack/concur.f',
  'fitpack/cualde.f',
  'fitpack/curev.f',
  'fitpack/curfit.f',
  'fitpack/dblint.f',
  'fitpack/evapol.f',
  'fitpack/fourco.f',
  'fitpack/fpader.f',
  'fitpack/fpadno.f',
  'fitpack/fpadpo.f',
  'fitpack/fpback.f',
  'fitpack/fpbacp.f',
  'fitpack/fpbfout.f',
  'fitpack/fpbisp.f',
  'fitpack/fpbspl.f',
  'fitpack/fpchec.f',
  'fitpack/fpched.f',
  'fitpack/fpchep.f',
  'fitpack/fpclos.f',
  'fitpack/fpcoco.f',
  'fitpack/fpcons.f',
  'fitpack/fpcosp.f',
  'fitpack/fpcsin.f',
  'fitpack/fpcurf.f',
  'fitpack/fpcuro.f',
  'fitpack/fpcyt1.f',
  'fitpack/fpcyt2.f',
  'fitpack/fpdeno.f',
  'fitpack/fpdisc.f',
  'fitpack/fpfrno.f',
  'fitpack/fpgivs.f',
  'fitpack/fpgrdi.f',
  'fitpack/fpgrpa.f',
  'fitpack/fpgrre.f',
  'fitpack/fpgrsp.f',
  'fitpack/fpinst.f',
  'fitpack/fpintb.f',
  'fitpack/fpknot.f',
  'fitpack/fpopdi.f',
  'fitpack/fpopsp.f',
  'fitpack/fporde.f',
  'fitpack/fppara.f',
  'fitpack/fppasu.f',
  'fitpack/fpperi.f',
  'fitpack/fppocu.f',
  'fitpack/fppogr.f',
  'fitpack/fppola.f',
  'fitpack/fprank.f',
  'fitpack/fprati.f',
  'fitpack/fpregr.f',
  'fitpack/fprota.f',
  'fitpack/fprppo.f',
  'fitpack/fprpsp.f',
  'fitpack/fpseno.f',
  'fitpack/fpspgr.f',
  'fitpack/fpsphe.f',
  'fitpack/fpsuev.f',
  'fitpack/fpsurf.f',
  'fitpack/fpsysy.f',
  'fitpack/fptrnp.f',
  'fitpack/fptrpe.f',
  'fitpack/insert.f',
  'fitpack/parcur.f',
  'fitpack/parder.f',
  'fitpack/pardeu.f',
  'fitpack/parsur.f',
  'fitpack/percur.f',
  'fitpack/pogrid.f',
  'fitpack/polar.f',
  'fitpack/profil.f',
  'fitpack/regrid.f',
  'fitpack/spalde.f',
  'fitpack/spgrid.f',
  'fitpack/sphere.f',
  'fitpack/splder.f',
  'fitpack/splev.f',
  'fitpack/splint.f',
  'fitpack/sproot.f',
  'fitpack/surev.f',
  'fitpack/surfit.f',
  'fitpack/pardtc.f'
]

# TODO: Add flags for 64 bit ints
fitpack_lib = static_library('fitpack_lib',
  fitpack_src,
  fortran_args: _fflag_Wno_maybe_uninitialized,
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

py3.extension_module('interpnd',
  spt_cython_gen.process('interpnd.pyx'),
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

py3.extension_module('_bspl',
  cython_gen.process('_bspl.pyx'),
  c_args: cython_c_args,
  include_directories: 'src/',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

# TODO: Add flags for 64 bit ints
py3.extension_module('_fitpack',
  ['src/_fitpackmodule.c'],
  link_with: [fitpack_lib],
  include_directories: 'src/',
  dependencies: np_dep,
  link_args: version_link_args,
  link_language: 'fortran',
  install: true,
  subdir: 'scipy/interpolate'
)

dfitpack_module = custom_target('dfitpack_module',
  output: ['dfitpack-f2pywrappers.f', 'dfitpackmodule.c'],
  input: 'src/fitpack.pyf',
  command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

# TODO: Add flags for 64 bit ints
py3.extension_module('dfitpack',
  dfitpack_module,
  c_args: [Wno_unused_variable],
  link_args: version_link_args,
  dependencies: [lapack_dep, fortranobject_dep],
  link_with: [fitpack_lib],
  override_options: ['b_lto=false'],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/interpolate'
)

if use_pythran
  py3.extension_module('_rbfinterp_pythran',
    pythran_gen.process('_rbfinterp_pythran.py'),
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/interpolate'
  )
else
  py3.install_sources(
    ['_rbfinterp_pythran.py'],
    subdir: 'scipy/interpolate'
  )
endif

py3.install_sources([
    '__init__.py',
    '_bsplines.py',
    '_cubic.py',
    '_fitpack_impl.py',
    '_fitpack_py.py',
    '_fitpack2.py',
    '_interpolate.py',
    '_ndgriddata.py',
    '_pade.py',
    '_polyint.py',
    '_rbf.py',
    '_rbfinterp.py',
    '_rgi.py',
    '_ndbspline.py',
    'fitpack.py',
    'fitpack2.py',
    'interpolate.py',
    'ndgriddata.py',
    'polyint.py',
    'rbf.py'
  ],
  subdir: 'scipy/interpolate'
)

subdir('tests')
