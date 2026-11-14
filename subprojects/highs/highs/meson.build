# -------------------------- Configuration
conf_data = configuration_data()
conf_data.set('FAST_BUILD',
              get_option('fast_build'))
conf_data.set('HIGHS_VERSION_MAJOR',
              meson.project_version().split('.')[0])
conf_data.set('HIGHS_VERSION_MINOR',
              meson.project_version().split('.')[1])
conf_data.set('HIGHS_VERSION_PATCH',
              meson.project_version().split('.')[2])
conf_data.set('ZLIB_FOUND',
              zlib_dep.found())
# Is the use of the following two lines the cause of the meson build
# failure?
conf_data.set('CUPDLP_CPU', true)
            #  cupdlp_cpu.found())
# 64 bit numbers
if host_machine.cpu_family() == 'x86_64'
    # Get user's option, if it's not provided, enable highsint64 by default on x86_64
    highsint64_opt = get_option('highsint64')
    conf_data.set('HIGHSINT64', highsint64_opt)
else
    conf_data.set('HIGHSINT64', false)
endif
# Date
# today_str = 'unknown'
# if (not meson.is_subproject())
#     python_getdate = '''
# import datetime
# import os
# import time

# (ig)
# Deprecate compilation date

# build_date = datetime.datetime.utcfromtimestamp(
#     int(os.environ.get('SOURCE_DATE_EPOCH', time.time()))
# )
# build_date_str = build_date.strftime('%Y-%m-%d')
# print(build_date_str)
#     '''
#     # Note: this is not guaranteed to work on Windows. This should be removed
#     #       (as proposed in gh-1735), and otherwise moved into a .py script with
#     #       a `#!/usr/bin/env python3` shebang, and `python3 -c` dropped.
#     today_cmd = run_command('python3', '-c',
#                             python_getdate,
#                             check: false)
#     today_str = today_cmd.stdout().strip()
# endif
# conf_data.set_quoted('HIGHS_COMPILATION_DATE',
#                      today_str)

conf_data.set_quoted('HIGHS_COMPILATION_DATE',
                     'deprecated')

# Conditional compiler options
_mm_pause_code = '''
#include <immintrin.h>
int main(){
  _mm_pause();
  return 0;
}
'''
_have_mm_pause = cppc.compiles(_mm_pause_code,
                              name: 'mm_pause check')
conf_data.set('HIGHS_HAVE_MM_PAUSE',
              _have_mm_pause)

if cppc.get_id() == 'msvc'
  _bitscan_rev_code = '''
  #include <intrin.h>
  #pragma intrinsic(_BitScanReverse)
  #pragma intrinsic(_BitScanReverse64)
  int main () {
   unsigned long x = 5;
   unsigned long y;
   _BitScanReverse(&y, x);
   _BitScanReverse64(&x, y);
   return 0;
   }
   '''
  _have_bitscan_reverse = cppc.compiles(_bitscan_rev_code,
                                       name: 'bitscan_rev check')
  conf_data.set('HIGHS_HAVE_BITSCAN_REVERSE',
                _have_bitscan_reverse)
else
  _builtin_clz_code = '''
  #include <cstdint>
  int main(){
    unsigned int x = 5;
    unsigned long long y = __builtin_clz(x);
    x = __builtin_clzll(y);
    return 0;
  }
  '''
  _have_builtin_clz = cppc.compiles(_builtin_clz_code,
                                       name: 'builtin_clz check')
  conf_data.set('HIGHS_HAVE_BUILTIN_CLZ',
                _have_builtin_clz)
endif

interim_conf = configure_file(
  input: 'HConfig.h.meson.in',
  output: 'HConfig.h.meson.interim',
  configuration: conf_data
)

highs_conf_file = vcs_tag(
  input: interim_conf,
  output: 'HConfig.h',
  fallback: 'unknown',
  replace_string: '_HIGHS_GITHASH_', # Remember to quote manually
)

_incdirs += include_directories('.')

# --------------------- Libraries

_cupdlp_srcs = [
    'pdlp/cupdlp/cupdlp_cs.c',
    'pdlp/cupdlp/cupdlp_linalg.c',
    'pdlp/cupdlp/cupdlp_proj.c',
    'pdlp/cupdlp/cupdlp_restart.c',
    'pdlp/cupdlp/cupdlp_scaling.c',
    'pdlp/cupdlp/cupdlp_solver.c',
    'pdlp/cupdlp/cupdlp_step.c',
    'pdlp/cupdlp/cupdlp_utils.c',
]

_basiclu_srcs = [
    'ipm/basiclu/basiclu_factorize.c',
    'ipm/basiclu/basiclu_get_factors.c',
    'ipm/basiclu/basiclu_initialize.c',
    'ipm/basiclu/basiclu_object.c',
    'ipm/basiclu/basiclu_solve_dense.c',
    'ipm/basiclu/basiclu_solve_for_update.c',
    'ipm/basiclu/basiclu_solve_sparse.c',
    'ipm/basiclu/basiclu_update.c',
    'ipm/basiclu/lu_build_factors.c',
    'ipm/basiclu/lu_condest.c',
    'ipm/basiclu/lu_dfs.c',
    'ipm/basiclu/lu_factorize_bump.c',
    'ipm/basiclu/lu_file.c',
    'ipm/basiclu/lu_garbage_perm.c',
    'ipm/basiclu/lu_initialize.c',
    'ipm/basiclu/lu_internal.c',
    'ipm/basiclu/lu_markowitz.c',
    'ipm/basiclu/lu_matrix_norm.c',
    'ipm/basiclu/lu_pivot.c',
    'ipm/basiclu/lu_residual_test.c',
    'ipm/basiclu/lu_setup_bump.c',
    'ipm/basiclu/lu_singletons.c',
    'ipm/basiclu/lu_solve_dense.c',
    'ipm/basiclu/lu_solve_for_update.c',
    'ipm/basiclu/lu_solve_sparse.c',
    'ipm/basiclu/lu_solve_symbolic.c',
    'ipm/basiclu/lu_solve_triangular.c',
    'ipm/basiclu/lu_update.c',
]

_ipx_srcs = [
    'ipm/ipx/basiclu_kernel.cc',
    'ipm/ipx/basiclu_wrapper.cc',
    'ipm/ipx/basis.cc',
    'ipm/ipx/conjugate_residuals.cc',
    'ipm/ipx/control.cc',
    'ipm/ipx/crossover.cc',
    'ipm/ipx/diagonal_precond.cc',
    'ipm/ipx/forrest_tomlin.cc',
    'ipm/ipx/guess_basis.cc',
    'ipm/ipx/indexed_vector.cc',
    'ipm/ipx/info.cc',
    'ipm/ipx/ipm.cc',
    'ipm/ipx/ipx_c.cc',
    'ipm/ipx/iterate.cc',
    'ipm/ipx/kkt_solver.cc',
    'ipm/ipx/kkt_solver_basis.cc',
    'ipm/ipx/kkt_solver_diag.cc',
    'ipm/ipx/linear_operator.cc',
    'ipm/ipx/lp_solver.cc',
    'ipm/ipx/lu_factorization.cc',
    'ipm/ipx/lu_update.cc',
    'ipm/ipx/maxvolume.cc',
    'ipm/ipx/model.cc',
    'ipm/ipx/normal_matrix.cc',
    'ipm/ipx/sparse_matrix.cc',
    'ipm/ipx/sparse_utils.cc',
    'ipm/ipx/splitted_normal_matrix.cc',
    'ipm/ipx/starting_basis.cc',
    'ipm/ipx/symbolic_invert.cc',
    'ipm/ipx/timer.cc',
    'ipm/ipx/utils.cc',
]

_hipo_srcs = [
    'ipm/hipo/ipm/IpmData.cpp',
    'ipm/hipo/ipm/FactorHighsSolver.cpp',
    'ipm/hipo/ipm/Control.cpp',
    'ipm/hipo/ipm/Iterate.cpp',
    'ipm/hipo/ipm/KktMatrix.cpp',
    'ipm/hipo/ipm/Model.cpp',
    'ipm/hipo/ipm/PreProcess.cpp',
    'ipm/hipo/ipm/Refine.cpp',
    'ipm/hipo/ipm/Solver.cpp',
    'ipm/hipo/ipm/UpLookingSolver.cpp',
    'ipm/hipo/factorhighs/Analyse.cpp',
    'ipm/hipo/factorhighs/CallAndTimeBlas.cpp',
    'ipm/hipo/factorhighs/CliqueStack.cpp',
    'ipm/hipo/factorhighs/DataCollector.cpp',
    'ipm/hipo/factorhighs/DenseFactHybrid.cpp',
    'ipm/hipo/factorhighs/DenseFactKernel.cpp',
    'ipm/hipo/factorhighs/DgemmParallel.cpp',
    'ipm/hipo/factorhighs/FactorHighs_c_api.cpp',
    'ipm/hipo/factorhighs/FactorHighs.cpp',
    'ipm/hipo/factorhighs/Factorise.cpp',
    'ipm/hipo/factorhighs/FormatHandler.cpp',
    'ipm/hipo/factorhighs/HybridHybridFormatHandler.cpp',
    'ipm/hipo/factorhighs/HybridSolveHandler.cpp',
    'ipm/hipo/factorhighs/KrylovMethodsIpm.cpp',
    'ipm/hipo/factorhighs/Numeric.cpp',
    'ipm/hipo/factorhighs/SolveHandler.cpp',
    'ipm/hipo/factorhighs/Swaps.cpp',
    'ipm/hipo/factorhighs/Symbolic.cpp',
    'ipm/hipo/auxiliary/Auxiliary.cpp',
    'ipm/hipo/auxiliary/KrylovMethods.cpp',
    'ipm/hipo/auxiliary/Logger.cpp',
    'ipm/hipo/auxiliary/VectorOperations.cpp',
]


_srcs = [
    'HighsExternalApi.cpp',
    'HighsExternalDeps.cpp',
    'interfaces/highs_c_api.cpp',
    'io/Filereader.cpp',
    'io/FilereaderLp.cpp',
    'io/FilereaderMps.cpp',
    'io/HMPSIO.cpp',
    'io/HMpsFF.cpp',
    'io/HighsIO.cpp',
    'io/LoadOptions.cpp',
    'io/LpReader.cpp',
    'ipm/IpxWrapper.cpp',
    'lp_data/Highs.cpp',
    'lp_data/HighsCallback.cpp',
    'lp_data/HighsDebug.cpp',
    'lp_data/HighsIis.cpp',
    'lp_data/HighsDeprecated.cpp',
    'lp_data/HighsInfo.cpp',
    'lp_data/HighsInfoDebug.cpp',
    'lp_data/HighsInterface.cpp',
    'lp_data/HighsLp.cpp',
    'lp_data/HighsLpUtils.cpp',
    'lp_data/HighsModelUtils.cpp',
    'lp_data/HighsOptions.cpp',
    'lp_data/HighsRanging.cpp',
    'lp_data/HighsRunData.cpp',
    'lp_data/HighsSolution.cpp',
    'lp_data/HighsSolutionDebug.cpp',
    'lp_data/HighsSolve.cpp',
    'lp_data/HighsStatus.cpp',
    'mip/HighsCliqueTable.cpp',
    'mip/HighsConflictPool.cpp',
    'mip/HighsCutGeneration.cpp',
    'mip/HighsCutPool.cpp',
    'mip/HighsDebugSol.cpp',
    'mip/HighsDomain.cpp',
    'mip/HighsDynamicRowMatrix.cpp',
    'mip/HighsFeasibilityJump.cpp',
    'mip/HighsGFkSolve.cpp',
    'mip/HighsImplications.cpp',
    'mip/HighsLpAggregator.cpp',
    'mip/HighsLpRelaxation.cpp',
    'mip/HighsMachineSchedSeparator.cpp',
    'mip/HighsMipSolver.cpp',
    'mip/HighsMipSolverData.cpp',
    'mip/HighsMipWorker.cpp',
    'mip/HighsModkSeparator.cpp',
    'mip/HighsNodeQueue.cpp',
    'mip/HighsObjectiveFunction.cpp',
    'mip/HighsPathSeparator.cpp',
    'mip/HighsPrimalHeuristics.cpp',
    'mip/HighsPseudocost.cpp',
    'mip/HighsRedcostFixing.cpp',
    'mip/HighsSearch.cpp',
    'mip/HighsSeparation.cpp',
    'mip/HighsSeparator.cpp',
    'mip/HighsTableauSeparator.cpp',
    'mip/HighsTransformedLp.cpp',
    'model/HighsHessian.cpp',
    'model/HighsHessianUtils.cpp',
    'model/HighsModel.cpp',
    'parallel/HighsTaskExecutor.cpp',
    'pdlp/CupdlpWrapper.cpp',
    'pdlp/HiPdlpWrapper.cpp',
    'pdlp/hipdlp/linalg.cc',
    'pdlp/hipdlp/logger.cc',
    'pdlp/hipdlp/pdhg.cc',
    'pdlp/hipdlp/restart.cc',
    'pdlp/hipdlp/scaling.cc',
    'presolve/HPresolve.cpp',
    'presolve/HPresolveAggregator.cpp',
    'presolve/HPresolveAnalysis.cpp',
    'presolve/HPresolveCliqueTable.cpp',
    'presolve/HPresolveColStuffing.cpp',
    'presolve/HPresolveDependent.cpp',
    'presolve/HPresolveDominatedCols.cpp',
    'presolve/HPresolveDualFixing.cpp',
    'presolve/HPresolveEnumeration.cpp',
    'presolve/HPresolveFourierMotzkin.cpp',
    'presolve/HPresolveImpliedInts.cpp',
    'presolve/HPresolveInitialSweep.cpp',
    'presolve/HPresolveParallelRowsAndCols.cpp',
    'presolve/HPresolveProbing.cpp',
    'presolve/HPresolveProbingUtil.cpp',
    'presolve/HPresolveSparsify.cpp',
    'presolve/HPresolveStrengthenIneq.cpp',
    'presolve/HPresolveUtils.cpp',
    'presolve/HPresolveVarBndAggregation.cpp',
    'presolve/HPresolveZeroCostSingleton.cpp',
    'presolve/HPresolveTest.cpp',
    'presolve/HighsPostsolveStack.cpp',
    'presolve/HighsSymmetry.cpp',
    'presolve/ICrash.cpp',
    'presolve/ICrashUtil.cpp',
    'presolve/ICrashX.cpp',
    'presolve/PresolveComponent.cpp',
    'qpsolver/a_asm.cpp',
    'qpsolver/a_quass.cpp',
    'qpsolver/basis.cpp',
    'qpsolver/QpAsmWrapper.cpp',
    'qpsolver/quass.cpp',
    'qpsolver/ratiotest.cpp',
    'qpsolver/scaling.cpp',
    'simplex/HEkk.cpp',
    'simplex/HEkkControl.cpp',
    'simplex/HEkkDebug.cpp',
    'simplex/HEkkDual.cpp',
    'simplex/HEkkDualMulti.cpp',
    'simplex/HEkkDualRHS.cpp',
    'simplex/HEkkDualRow.cpp',
    'simplex/HEkkInterface.cpp',
    'simplex/HEkkPrimal.cpp',
    'simplex/HSimplex.cpp',
    'simplex/HSimplexDebug.cpp',
    'simplex/HSimplexNla.cpp',
    'simplex/HSimplexNlaDebug.cpp',
    'simplex/HSimplexNlaFreeze.cpp',
    'simplex/HSimplexNlaProductForm.cpp',
    'simplex/HSimplexReport.cpp',
    'simplex/HighsSimplexAnalysis.cpp',
    'util/HFactor.cpp',
    'util/HFactorDebug.cpp',
    'util/HFactorExtend.cpp',
    'util/HFactorRefactor.cpp',
    'util/HFactorUtils.cpp',
    'util/HSet.cpp',
    'util/HVectorBase.cpp',
    'util/HighsDynamicLibrary.cpp',
    'util/HighsHash.cpp',
    'util/HighsLinearSumBounds.cpp',
    'util/HighsMatrixPic.cpp',
    'util/HighsMatrixUtils.cpp',
    'util/HighsSort.cpp',
    'util/HighsSparseMatrix.cpp',
    'util/HighsUtils.cpp',
    'util/stringutil.cpp',
]

highslib_srcs = [
  files(_srcs),
  files(_cupdlp_srcs),
  files(_basiclu_srcs),
  files(_ipx_srcs),
  files(_hipo_srcs),
]


symbol_visibility = 'default'
if get_option('default_library') == 'static'
  # Ensure that if we link a static library into a shared library,
  # private symbols don't get re-exported.
  symbol_visibility = 'inlineshidden'
endif

do_install = get_option('do_install')

if not get_option('with_pybind11')
highslib = library('highs',
                   highslib_srcs + highs_conf_file,
                   dependencies: _deps,
                   cpp_args: _highs_cpp_suppressions + _args,
                   c_args: _highs_c_suppressions + _args,
                   include_directories: _incdirs,
                   gnu_symbol_visibility: symbol_visibility,
                   pic: true,
                   install: do_install)

highs_dep = declare_dependency(link_with: highslib,
                               dependencies: _deps,
                               include_directories: _incdirs,
                               )
endif

# --------------- Interfaces

# Fortran
if get_option('with_fortran')
  add_languages('fortran', required: true)
  _fsrcs = [
    'interfaces/highs_fortran_api.f90'
  ]
  library('FortranHighs',
          _fsrcs,
          dependencies: _deps,
          cpp_args: _args,
          link_with: highslib,
          include_directories: _incdirs,
          pic: true,
          install: do_install)
endif

# C#
if get_option('with_csharp')
  add_languages('cs', required: true)
  _cs_srcs = [
    'interfaces/highs_csharp_api.cs'
  ]
  library('HighsCsharp',
          _cs_srcs,
          dependencies: _deps,
          cpp_args: _args,
          link_with: highslib,
          include_directories: _incdirs,
          pic: true,
          install: do_install)
endif

highspy_cpp = files([
  '../highspy/highs_bindings.cpp'
])

if get_option('with_pybind11')
  py = import('python').find_installation(pure: false)

  highslib = library('highs',
                     highslib_srcs + highs_conf_file,
                     dependencies: _deps,
                     cpp_args: _highs_cpp_suppressions + _args,
                     c_args: _highs_c_suppressions + _args,
                     include_directories: _incdirs,
                     gnu_symbol_visibility: symbol_visibility,
                     pic: true,
                     install: do_install,
                     install_dir: py.get_install_dir() / 'highspy')

  highs_dep = declare_dependency(link_with: highslib,
                                 dependencies: _deps,
                                 include_directories: _incdirs,
                                )
  _deps += highs_dep
  pyb11_dep = [
    # py_dep is auto-added for Python >= 3.9, so it can be dropped here when
    # that is the minimum supported Python version
    py.dependency(),
    dependency('pybind11')
  ]
  _deps += pyb11_dep

  highspy_cpp = files([
    'highs_bindings.cpp'
  ])

  py.extension_module(
    '_core',
    sources : highspy_cpp,
    dependencies: [pyb11_dep, highs_dep],
    cpp_args: _args,
    install: do_install,
    subdir: 'highspy',
    include_directories: _incdirs,
  )

 highs_py_srcs = files(
      'highspy/__init__.py',
      'highspy/highs.py',)

  py.install_sources(highs_py_srcs,
    pure: false,
    subdir: 'highspy'
  )
endif
