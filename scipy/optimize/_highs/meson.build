highs_define_macros = [
  '-DCMAKE_BUILD_TYPE="RELEASE"',
  '-DHiGHSRELEASE',
  '-DIPX_ON="ON"',
  '-DHIGHS_GITHASH="n/a"',
  '-DHIGHS_COMPILATION_DATE="2021-07-09"',  # cannot generate dynamically
  '-DHIGHS_VERSION_MAJOR=1',  # don't care about this, look at CMakelists.txt
  '-DHIGHS_VERSION_MINOR=0',
  '-DHIGHS_VERSION_PATCH=0',
  '-DHIGHS_DIR=' + meson.current_source_dir()
]

basiclu_lib = static_library('basiclu',
  [
    'src/ipm/basiclu/src/basiclu_factorize.c',
    'src/ipm/basiclu/src/basiclu_get_factors.c',
    'src/ipm/basiclu/src/basiclu_initialize.c',
    'src/ipm/basiclu/src/basiclu_object.c',
    'src/ipm/basiclu/src/basiclu_solve_dense.c',
    'src/ipm/basiclu/src/basiclu_solve_for_update.c',
    'src/ipm/basiclu/src/basiclu_solve_sparse.c',
    'src/ipm/basiclu/src/basiclu_update.c',
    'src/ipm/basiclu/src/lu_build_factors.c',
    'src/ipm/basiclu/src/lu_condest.c',
    'src/ipm/basiclu/src/lu_dfs.c',
    'src/ipm/basiclu/src/lu_factorize_bump.c',
    'src/ipm/basiclu/src/lu_file.c',
    'src/ipm/basiclu/src/lu_garbage_perm.c',
    'src/ipm/basiclu/src/lu_initialize.c',
    'src/ipm/basiclu/src/lu_internal.c',
    'src/ipm/basiclu/src/lu_markowitz.c',
    'src/ipm/basiclu/src/lu_matrix_norm.c',
    'src/ipm/basiclu/src/lu_pivot.c',
    'src/ipm/basiclu/src/lu_residual_test.c',
    'src/ipm/basiclu/src/lu_setup_bump.c',
    'src/ipm/basiclu/src/lu_singletons.c',
    'src/ipm/basiclu/src/lu_solve_dense.c',
    'src/ipm/basiclu/src/lu_solve_for_update.c',
    'src/ipm/basiclu/src/lu_solve_sparse.c',
    'src/ipm/basiclu/src/lu_solve_symbolic.c',
    'src/ipm/basiclu/src/lu_solve_triangular.c',
    'src/ipm/basiclu/src/lu_update.c'
  ],
  include_directories: ['src', 'src/ipm/basiclu/include'],
  c_args: [
    '-Wno-unused-variable', highs_define_macros
  ]
)

# highs_wrapper
ipx_sources = [
  'src/ipm/ipx/src/basiclu_kernel.cc',
  'src/ipm/ipx/src/basiclu_wrapper.cc',
  'src/ipm/ipx/src/basis.cc',
  'src/ipm/ipx/src/conjugate_residuals.cc',
  'src/ipm/ipx/src/control.cc',
  'src/ipm/ipx/src/crossover.cc',
  'src/ipm/ipx/src/diagonal_precond.cc',
  'src/ipm/ipx/src/forrest_tomlin.cc',
  'src/ipm/ipx/src/guess_basis.cc',
  'src/ipm/ipx/src/indexed_vector.cc',
  'src/ipm/ipx/src/info.cc',
  'src/ipm/ipx/src/ipm.cc',
  'src/ipm/ipx/src/ipx_c.cc',
  'src/ipm/ipx/src/iterate.cc',
  'src/ipm/ipx/src/kkt_solver.cc',
  'src/ipm/ipx/src/kkt_solver_basis.cc',
  'src/ipm/ipx/src/kkt_solver_diag.cc',
  'src/ipm/ipx/src/linear_operator.cc',
  'src/ipm/ipx/src/lp_solver.cc',
  'src/ipm/ipx/src/lu_factorization.cc',
  'src/ipm/ipx/src/lu_update.cc',
  'src/ipm/ipx/src/maxvolume.cc',
  'src/ipm/ipx/src/model.cc',
  'src/ipm/ipx/src/normal_matrix.cc',
  'src/ipm/ipx/src/sparse_matrix.cc',
  'src/ipm/ipx/src/sparse_utils.cc',
  'src/ipm/ipx/src/splitted_normal_matrix.cc',
  'src/ipm/ipx/src/starting_basis.cc',
  'src/ipm/ipx/src/symbolic_invert.cc',
  'src/ipm/ipx/src/timer.cc',
  'src/ipm/ipx/src/utils.cc'
]

highs_sources = [
  'external/filereaderlp/reader.cpp',
  'src/interfaces/highs_c_api.cpp',
  'src/io/Filereader.cpp',
  'src/io/FilereaderEms.cpp',
  'src/io/FilereaderLp.cpp',
  'src/io/FilereaderMps.cpp',
  'src/io/HMPSIO.cpp',
  'src/io/HMpsFF.cpp',
  'src/io/HighsIO.cpp',
  'src/io/LoadOptions.cpp',
  'src/lp_data/Highs.cpp',
  'src/lp_data/HighsDebug.cpp',
  'src/lp_data/HighsInfo.cpp',
  'src/lp_data/HighsLp.cpp',
  'src/lp_data/HighsLpUtils.cpp',
  'src/lp_data/HighsModelUtils.cpp',
  'src/lp_data/HighsOptions.cpp',
  'src/lp_data/HighsRanging.cpp',
  'src/lp_data/HighsSolution.cpp',
  'src/lp_data/HighsSolutionDebug.cpp',
  'src/lp_data/HighsSolve.cpp',
  'src/lp_data/HighsStatus.cpp',
  'src/presolve/HPreData.cpp',
  'src/presolve/Presolve.cpp',
  'src/presolve/PresolveAnalysis.cpp',
  'src/presolve/PresolveComponent.cpp',
  'src/presolve/PresolveUtils.cpp',
  'src/simplex/HCrash.cpp',
  'src/simplex/HDual.cpp',
  'src/simplex/HDualMulti.cpp',
  'src/simplex/HDualRHS.cpp',
  'src/simplex/HDualRow.cpp',
  'src/simplex/HFactor.cpp',
  'src/simplex/HFactorDebug.cpp',
  'src/simplex/HMatrix.cpp',
  'src/simplex/HPrimal.cpp',
  'src/simplex/HQPrimal.cpp',
  'src/simplex/HSimplex.cpp',
  'src/simplex/HSimplexDebug.cpp',
  'src/simplex/HSimplexReport.cpp',
  'src/simplex/HVector.cpp',
  'src/simplex/HighsSimplexAnalysis.cpp',
  'src/simplex/HighsSimplexInterface.cpp',
  'src/test/DevKkt.cpp',
  'src/test/KktCh2.cpp',
  'src/util/HighsMatrixPic.cpp',
  'src/util/HighsSort.cpp',
  'src/util/HighsUtils.cpp',
  'src/util/stringutil.cpp'
]

if cc.has_argument('-Wno-unused-but-set-variable')
  Wno_unused_but_set = ['-Wno-unused-but-set-variable']
else
  Wno_unused_but_set = []
endif

_highs_wrapper = py3.extension_module('_highs_wrapper',
  ['cython/src/_highs_wrapper.pyx', highs_sources, ipx_sources],
  include_directories: [
    inc_np,
    # highs_wrapper
    'cython/src/',
    'src/',
    'src/lp_data/',
    # highs
    'src/',
    'src/io/',
    'src/ipm/ipx/include/',
    # IPX
    'src/ipm/ipx/include/',
    'src/ipm/basiclu/include/',
  ],
  link_with: [basiclu_lib],
  dependencies: [py3_dep],
  cpp_args: [
    '-Wno-unused-variable',
    Wno_unused_but_set,
    highs_define_macros,
    '-UOPENMP',
    '-UEXT_PRESOLVE',
    '-USCIP_DEV',
    '-UHiGHSDEV',
    '-UOSI_FOUND',
    numpy_nodepr_api,
  ],
  override_options: ['cython_language=cpp'],
  install: true,
  subdir: 'scipy/optimize/_highs'
)

_highs_constants = py3.extension_module('_highs_constants',
  'cython/src/_highs_constants.pyx',
  include_directories: [
    'cython/src/',
    'src/',
    'src/io/',
    'src/lp_data/',
    'src/simplex/'
  ],
  dependencies: [py3_dep],
  override_options: ['cython_language=cpp'],
  install: true,
  subdir: 'scipy/optimize/_highs'
)

py3.install_sources([
    'cython/src/HConst.pxd',
    'cython/src/Highs.pxd',
    'cython/src/HighsIO.pxd',
    'cython/src/HighsInfo.pxd',
    'cython/src/HighsLp.pxd',
    'cython/src/HighsLpUtils.pxd',
    'cython/src/HighsModelUtils.pxd',
    'cython/src/HighsOptions.pxd',
    'cython/src/HighsRuntimeOptions.pxd',
    'cython/src/HighsStatus.pxd',
    'cython/src/SimplexConst.pxd',
    'cython/src/highs_c_api.pxd'
  ],
  pure: false,
  subdir: 'scipy/optimize/_highs/src/cython'
)

py3.install_sources(
  ['__init__.py'],
  pure: false,
  subdir: 'scipy/optimize/_highs'
)
