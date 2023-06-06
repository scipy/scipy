"""
setup.py for HiGHS scipy interface
"""

from os.path import join

from scipy._lib._highs_utils import _highs_dir


_highs_flags = [
    '-Wno-class-memaccess',
    '-Wno-format-truncation',
    '-Wno-non-virtual-dtor',
    '-Wno-sign-compare',
    '-Wno-switch',
    '-Wno-unused-but-set-variable',
    '-Wno-unused-variable',
]


def highs_clib_hook_cxx_flags(build_clib, build_info):
    from scipy._build_utils.compiler_helper import (set_cxx_flags_clib_hook,
                                                    try_add_flag)
    cc = build_clib.compiler
    if cc.compiler_type != 'msvc':
        if 'extra_compiler_args' not in build_info:
            build_info['extra_compiler_args'] = []
        for flag in _highs_flags:
            try_add_flag(build_info['extra_compiler_args'], cc, flag)
    set_cxx_flags_clib_hook(build_clib, build_info)


def pre_build_hook(build_ext, ext):
    from scipy._build_utils.compiler_helper import (set_cxx_flags_hook,
                                                    try_add_flag)
    cc = build_ext._cxx_compiler
    args = ext.extra_compile_args
    set_cxx_flags_hook(build_ext, ext)
    if cc.compiler_type != 'msvc':
        for flag in _highs_flags:
            try_add_flag(args, cc, flag)


def basiclu_pre_build_hook(build_clib, build_info):
    from scipy._build_utils.compiler_helper import get_c_std_flag
    c_flag = get_c_std_flag(build_clib.compiler)
    if c_flag is not None:
        if 'extra_compiler_args' not in build_info:
            build_info['extra_compiler_args'] = []
        build_info['extra_compiler_args'].append(c_flag)


def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    from scipy._build_utils import numpy_nodepr_api
    import pybind11
    config = Configuration('_highs', parent_package, top_path)

    # Here are the pound defines that HConfig.h would usually provide;
    # We provide an empty HConfig.h file and do the defs and undefs
    # here;  should be consistent with meson.build in this directory
    DEFINE_MACROS = [
        ('CMAKE_BUILD_TYPE', '"RELEASE"'),
        ('FAST_BUILD', 'ON'),
        ('HIGHS_GITHASH', '"n/a"'),
        ('HIGHS_COMPILATION_DATE', '"0000-00-00"'),
        ('HIGHS_VERSION_MAJOR', "1"),
        ('HIGHS_VERSION_MINOR', "5"),
        ('HIGHS_VERSION_PATCH', "3"),
        ('HIGHS_DIR', '"' + str(_highs_dir().resolve()) + '"'),
    ]
    UNDEF_MACROS = [
        'OPENMP',
        'EXT_PRESOLVE',
        'SCIP_DEV',
        'HiGHSDEV',
        'OSI_FOUND',
        'NDEBUG',
    ]

    # Compile BASICLU as a static library to appease clang:
    # (won't allow -std=c++11/14 option for C sources)
    basiclu_sources = [
        join('src', 'ipm', 'basiclu', 'basiclu_factorize.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_solve_dense.c'),
        join('src', 'ipm', 'basiclu', 'lu_build_factors.c'),
        join('src', 'ipm', 'basiclu', 'lu_factorize_bump.c'),
        join('src', 'ipm', 'basiclu', 'lu_initialize.c'),
        join('src', 'ipm', 'basiclu', 'lu_markowitz.c'),
        join('src', 'ipm', 'basiclu', 'lu_setup_bump.c'),
        join('src', 'ipm', 'basiclu', 'lu_solve_sparse.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_get_factors.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_solve_for_update.c'),
        join('src', 'ipm', 'basiclu', 'lu_condest.c'),
        join('src', 'ipm', 'basiclu', 'lu_file.c'),
        join('src', 'ipm', 'basiclu', 'lu_internal.c'),
        join('src', 'ipm', 'basiclu', 'lu_matrix_norm.c'),
        join('src', 'ipm', 'basiclu', 'lu_singletons.c'),
        join('src', 'ipm', 'basiclu', 'lu_solve_symbolic.c'),
        join('src', 'ipm', 'basiclu', 'lu_update.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_initialize.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_solve_sparse.c'),
        join('src', 'ipm', 'basiclu', 'lu_pivot.c'),
        join('src', 'ipm', 'basiclu', 'lu_solve_dense.c'),
        join('src', 'ipm', 'basiclu', 'lu_solve_triangular.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_object.c'),
        join('src', 'ipm', 'basiclu', 'basiclu_update.c'),
        join('src', 'ipm', 'basiclu', 'lu_dfs.c'),
        join('src', 'ipm', 'basiclu', 'lu_garbage_perm.c'),
        join('src', 'ipm', 'basiclu', 'lu_residual_test.c'),
        join('src', 'ipm', 'basiclu', 'lu_solve_for_update.c')
    ]

    ipx_sources = [
        join('src', 'ipm', 'ipx', 'basiclu_kernel.cc'),
        join('src', 'ipm', 'ipx', 'basiclu_wrapper.cc'),
        join('src', 'ipm', 'ipx', 'basis.cc'),
        join('src', 'ipm', 'ipx', 'conjugate_residuals.cc'),
        join('src', 'ipm', 'ipx', 'control.cc'),
        join('src', 'ipm', 'ipx', 'crossover.cc'),
        join('src', 'ipm', 'ipx', 'diagonal_precond.cc'),
        join('src', 'ipm', 'ipx', 'forrest_tomlin.cc'),
        join('src', 'ipm', 'ipx', 'guess_basis.cc'),
        join('src', 'ipm', 'ipx', 'indexed_vector.cc'),
        join('src', 'ipm', 'ipx', 'info.cc'),
        join('src', 'ipm', 'ipx', 'ipm.cc'),
        join('src', 'ipm', 'ipx', 'ipx_c.cc'),
        join('src', 'ipm', 'ipx', 'iterate.cc'),
        join('src', 'ipm', 'ipx', 'kkt_solver.cc'),
        join('src', 'ipm', 'ipx', 'kkt_solver_basis.cc'),
        join('src', 'ipm', 'ipx', 'kkt_solver_diag.cc'),
        join('src', 'ipm', 'ipx', 'linear_operator.cc'),
        join('src', 'ipm', 'ipx', 'lp_solver.cc'),
        join('src', 'ipm', 'ipx', 'lu_factorization.cc'),
        join('src', 'ipm', 'ipx', 'lu_update.cc'),
        join('src', 'ipm', 'ipx', 'maxvolume.cc'),
        join('src', 'ipm', 'ipx', 'model.cc'),
        join('src', 'ipm', 'ipx', 'normal_matrix.cc'),
        join('src', 'ipm', 'ipx', 'sparse_matrix.cc'),
        join('src', 'ipm', 'ipx', 'sparse_utils.cc'),
        join('src', 'ipm', 'ipx', 'splitted_normal_matrix.cc'),
        join('src', 'ipm', 'ipx', 'starting_basis.cc'),
        join('src', 'ipm', 'ipx', 'symbolic_invert.cc'),
        join('src', 'ipm', 'ipx', 'timer.cc'),
        join('src', 'ipm', 'ipx', 'utils.cc')
    ]
    highs_sources = [
        join('extern', 'filereaderlp', 'reader.cpp'),
        join('src', 'io', 'Filereader.cpp'),
        join('src', 'io', 'FilereaderLp.cpp'),
        join('src', 'io', 'FilereaderEms.cpp'),
        join('src', 'io', 'FilereaderMps.cpp'),
        join('src', 'io', 'HighsIO.cpp'),
        join('src', 'io', 'HMPSIO.cpp'),
        join('src', 'io', 'HMpsFF.cpp'),
        join('src', 'io', 'LoadOptions.cpp'),
        join('src', 'lp_data', 'Highs.cpp'),
        join('src', 'lp_data', 'HighsDebug.cpp'),
        join('src', 'lp_data', 'HighsDeprecated.cpp'),
        join('src', 'lp_data', 'HighsInfo.cpp'),
        join('src', 'lp_data', 'HighsInfoDebug.cpp'),
        join('src', 'lp_data', 'HighsInterface.cpp'),
        join('src', 'lp_data', 'HighsLp.cpp'),
        join('src', 'lp_data', 'HighsLpUtils.cpp'),
        join('src', 'lp_data', 'HighsModelUtils.cpp'),
        join('src', 'lp_data', 'HighsRanging.cpp'),
        join('src', 'lp_data', 'HighsSolution.cpp'),
        join('src', 'lp_data', 'HighsSolutionDebug.cpp'),
        join('src', 'lp_data', 'HighsSolve.cpp'),
        join('src', 'lp_data', 'HighsStatus.cpp'),
        join('src', 'lp_data', 'HighsOptions.cpp'),
        join('src', 'presolve', 'ICrash.cpp'),
        join('src', 'presolve', 'ICrashUtil.cpp'),
        join('src', 'presolve', 'ICrashX.cpp'),
        join('src', 'mip', 'HighsMipSolver.cpp'),
        join('src', 'mip', 'HighsMipSolverData.cpp'),
        join('src', 'mip', 'HighsDomain.cpp'),
        join('src', 'mip', 'HighsDynamicRowMatrix.cpp'),
        join('src', 'mip', 'HighsLpRelaxation.cpp'),
        join('src', 'mip', 'HighsSeparation.cpp'),
        join('src', 'mip', 'HighsSeparator.cpp'),
        join('src', 'mip', 'HighsTableauSeparator.cpp'),
        join('src', 'mip', 'HighsModkSeparator.cpp'),
        join('src', 'mip', 'HighsPathSeparator.cpp'),
        join('src', 'mip', 'HighsCutGeneration.cpp'),
        join('src', 'mip', 'HighsSearch.cpp'),
        join('src', 'mip', 'HighsConflictPool.cpp'),
        join('src', 'mip', 'HighsCutPool.cpp'),
        join('src', 'mip', 'HighsCliqueTable.cpp'),
        join('src', 'mip', 'HighsGFkSolve.cpp'),
        join('src', 'mip', 'HighsTransformedLp.cpp'),
        join('src', 'mip', 'HighsLpAggregator.cpp'),
        join('src', 'mip', 'HighsDebugSol.cpp'),
        join('src', 'mip', 'HighsImplications.cpp'),
        join('src', 'mip', 'HighsPrimalHeuristics.cpp'),
        join('src', 'mip', 'HighsPseudocost.cpp'),
        join('src', 'mip', 'HighsNodeQueue.cpp'),
        join('src', 'mip', 'HighsObjectiveFunction.cpp'),
        join('src', 'mip', 'HighsRedcostFixing.cpp'),
        join('src', 'model', 'HighsHessian.cpp'),
        join('src', 'model', 'HighsHessianUtils.cpp'),
        join('src', 'model', 'HighsModel.cpp'),
        join('src', 'parallel', 'HighsTaskExecutor.cpp'),
        join('src', 'presolve', 'HighsPostsolveStack.cpp'),
        join('src', 'presolve', 'HighsSymmetry.cpp'),
        join('src', 'presolve', 'HPresolve.cpp'),
        join('src', 'presolve', 'HPresolveAnalysis.cpp'),
        join('src', 'presolve', 'PresolveComponent.cpp'),
        join('src', 'qpsolver', 'basis.cpp'),
        join('src', 'qpsolver', 'quass.cpp'),
        join('src', 'qpsolver', 'ratiotest.cpp'),
        join('src', 'qpsolver', 'scaling.cpp'),
        join('src', 'qpsolver', 'perturbation.cpp'),
        join('src', 'simplex', 'HEkk.cpp'),
        join('src', 'simplex', 'HEkkControl.cpp'),
        join('src', 'simplex', 'HEkkDebug.cpp'),
        join('src', 'simplex', 'HEkkPrimal.cpp'),
        join('src', 'simplex', 'HEkkDual.cpp'),
        join('src', 'simplex', 'HEkkDualRHS.cpp'),
        join('src', 'simplex', 'HEkkDualRow.cpp'),
        join('src', 'simplex', 'HEkkDualMulti.cpp'),
        join('src', 'simplex', 'HEkkInterface.cpp'),
        join('src', 'simplex', 'HighsSimplexAnalysis.cpp'),
        join('src', 'simplex', 'HSimplex.cpp'),
        join('src', 'simplex', 'HSimplexDebug.cpp'),
        join('src', 'simplex', 'HSimplexNla.cpp'),
        join('src', 'simplex', 'HSimplexNlaDebug.cpp'),
        join('src', 'simplex', 'HSimplexNlaFreeze.cpp'),
        join('src', 'simplex', 'HSimplexNlaProductForm.cpp'),
        join('src', 'simplex', 'HSimplexReport.cpp'),
        join('src', 'test', 'KktCh2.cpp'),
        join('src', 'test', 'DevKkt.cpp'),
        join('src', 'util', 'HFactor.cpp'),
        join('src', 'util', 'HFactorDebug.cpp'),
        join('src', 'util', 'HFactorExtend.cpp'),
        join('src', 'util', 'HFactorRefactor.cpp'),
        join('src', 'util', 'HFactorUtils.cpp'),
        join('src', 'util', 'HighsHash.cpp'),
        join('src', 'util', 'HighsLinearSumBounds.cpp'),
        join('src', 'util', 'HighsMatrixPic.cpp'),
        join('src', 'util', 'HighsMatrixUtils.cpp'),
        join('src', 'util', 'HighsSort.cpp'),
        join('src', 'util', 'HighsSparseMatrix.cpp'),
        join('src', 'util', 'HighsUtils.cpp'),
        join('src', 'util', 'HSet.cpp'),
        join('src', 'util', 'HVectorBase.cpp'),
        join('src', 'util', 'stringutil.cpp'),
        join('src', 'interfaces', 'highs_c_api.cpp'),
        join('src', 'ipm', 'IpxWrapper.cpp')
    ]

    basiclu_sources = [str(_highs_dir() / s) for s in basiclu_sources]
    ipx_sources = [str(_highs_dir() / s) for s in ipx_sources]
    highs_sources = [str(_highs_dir() / s) for s in highs_sources]

    # Compile BASICLU as a static library to appease clang:
    # (won't allow -std=c++11/14 option for C sources)
    root = str(_highs_dir().resolve())
    config.add_library(
        'basiclu',
        sources=basiclu_sources,
        include_dirs=[
            'src/',
            join(root, 'src'),
            join(root, 'src', 'util'),
            join(root, 'extern'),
            join(root, 'src', 'ipm', 'basiclu'),
        ],
        language='c',
        macros=DEFINE_MACROS,
        _pre_build_hook=basiclu_pre_build_hook,
    )

    config.add_library(
        'ipx',
        sources=ipx_sources,
        libraries=['basiclu'],
        include_dirs=[
            'src',
            join(root, 'extern'),
            join(root, 'src'),
            join(root, 'src', 'ipm', 'ipx'),
        ],
        language='c++',
        macros=DEFINE_MACROS,
        _pre_build_hook=highs_clib_hook_cxx_flags,
    )

    config.add_library(
        'highs',
        sources=highs_sources,
        libraries=['ipx', 'basiclu'],
        include_dirs=[
            'src',
            join(root, 'extern'),
            join(root, 'src'),
            join(root, 'src', 'io'),
            join(root, 'src', 'ipm', 'ipx'),
            join(root, 'src', 'lp_data'),
            join(root, 'src', 'util'),
        ],
        language='c++',
        macros=DEFINE_MACROS,
        _pre_build_hook=highs_clib_hook_cxx_flags,
    )

    highs_bindings_pybind_includes = [
        pybind11.get_include(True),
        pybind11.get_include(False),
        get_numpy_include_dirs(),
    ]
    ext = config.add_extension(
        'highs_bindings',
        sources=[join(root, 'highspy', 'highs_bindings.cpp')],
        define_macros=DEFINE_MACROS + numpy_nodepr_api['define_macros'],
        undef_macros=UNDEF_MACROS,
        include_dirs=[
            'src',
            join(root, 'extern'),
            join(root, 'src'),
            join(root, 'src', 'io'),
            join(root, 'src', 'ipm', 'ipx'),
            join(root, 'src', 'lp_data'),
            join(root, 'src', 'util'),
        ] + highs_bindings_pybind_includes,
        libraries=['highs', 'ipx', 'basiclu'],
        language='c++',
    )
    ext._pre_build_hook = pre_build_hook

    ext = config.add_extension(
        '_highs_options',
        sources=['_highs_options.cpp'],
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        include_dirs=[
            'src',
            join(root, 'src'),
            join(root, 'src', 'lp_data')
        ],
        libraries=['highs'],
        language='c++',
    )
    ext._pre_build_hook = pre_build_hook

    ext = config.add_extension(
        '_highs_constants',
        sources=['_highs_constants.cpp'],
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        include_dirs=[
            'src',
            join(root, 'src'),
            join(root, 'src', 'lp_data'),
            join(root, 'src', 'simplex')
        ],
        libraries=['highs'],
        language='c++',
    )
    ext._pre_build_hook = pre_build_hook

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
