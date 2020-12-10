# Define library.
# Outdated CMake approach: update in progress
if (NOT FAST_BUILD)
include_directories(ipm/ipx/src)
include_directories(ipm/ipx/include)
include_directories(ipm/basiclu/src)
include_directories(ipm/basiclu/include)

if (IPX_ON)
set(basiclu_sources
    ipm/basiclu/src/basiclu_factorize.c
    ipm/basiclu/src/basiclu_solve_dense.c
    ipm/basiclu/src/lu_build_factors.c
    ipm/basiclu/src/lu_factorize_bump.c
    ipm/basiclu/src/lu_initialize.c
    ipm/basiclu/src/lu_markowitz.c
    ipm/basiclu/src/lu_setup_bump.c
    ipm/basiclu/src/lu_solve_sparse.c
    ipm/basiclu/src/basiclu_get_factors.c
    ipm/basiclu/src/basiclu_solve_for_update.c
    ipm/basiclu/src/lu_condest.c
    ipm/basiclu/src/lu_file.c
    ipm/basiclu/src/lu_internal.c
    ipm/basiclu/src/lu_matrix_norm.c
    ipm/basiclu/src/lu_singletons.c
    ipm/basiclu/src/lu_solve_symbolic.c
    ipm/basiclu/src/lu_update.c
    ipm/basiclu/src/basiclu_initialize.c
    ipm/basiclu/src/basiclu_solve_sparse.c
    ipm/basiclu/src/lu_pivot.c
    ipm/basiclu/src/lu_solve_dense.c
    ipm/basiclu/src/lu_solve_triangular.c
    ipm/basiclu/src/basiclu_object.c
    ipm/basiclu/src/basiclu_update.c
    ipm/basiclu/src/lu_dfs.c
    ipm/basiclu/src/lu_garbage_perm.c
    ipm/basiclu/src/lu_residual_test.c
    ipm/basiclu/src/lu_solve_for_update.c)

set(ipx_sources
    ipm/ipx/src/basiclu_kernel.cc
    ipm/ipx/src/basiclu_wrapper.cc
    ipm/ipx/src/basis.cc
    ipm/ipx/src/conjugate_residuals.cc
    ipm/ipx/src/control.cc
    ipm/ipx/src/crossover.cc
    ipm/ipx/src/diagonal_precond.cc
    ipm/ipx/src/forrest_tomlin.cc
    ipm/ipx/src/guess_basis.cc
    ipm/ipx/src/indexed_vector.cc
    ipm/ipx/src/info.cc
    ipm/ipx/src/ipm.cc
    ipm/ipx/src/ipx_c.cc
    ipm/ipx/src/iterate.cc
    ipm/ipx/src/kkt_solver.cc
    ipm/ipx/src/kkt_solver_basis.cc
    ipm/ipx/src/kkt_solver_diag.cc
    ipm/ipx/src/linear_operator.cc
    ipm/ipx/src/lp_solver.cc
    ipm/ipx/src/lu_factorization.cc
    ipm/ipx/src/lu_update.cc
    ipm/ipx/src/maxvolume.cc
    ipm/ipx/src/model.cc
    ipm/ipx/src/normal_matrix.cc
    ipm/ipx/src/sparse_matrix.cc
    ipm/ipx/src/sparse_utils.cc
    ipm/ipx/src/splitted_normal_matrix.cc
    ipm/ipx/src/starting_basis.cc
    ipm/ipx/src/symbolic_invert.cc
    ipm/ipx/src/timer.cc
    ipm/ipx/src/utils.cc)
endif()

set(sources
    ../external/filereaderlp/reader.cpp
    io/Filereader.cpp
    io/FilereaderLp.cpp
    io/FilereaderEms.cpp
    io/FilereaderMps.cpp
    io/HighsIO.cpp
    io/HMPSIO.cpp
    io/HMpsFF.cpp
    io/LoadOptions.cpp
    lp_data/Highs.cpp
    lp_data/HighsDebug.cpp
    lp_data/HighsInfo.cpp
    lp_data/HighsLp.cpp
    lp_data/HighsLpUtils.cpp
    lp_data/HighsModelUtils.cpp
    lp_data/HighsRanging.cpp
    lp_data/HighsSolution.cpp
    lp_data/HighsSolutionDebug.cpp
    lp_data/HighsSolve.cpp
    lp_data/HighsStatus.cpp
    lp_data/HighsOptions.cpp
    mip/HighsMipSolver.cpp
    mip/SolveMip.cpp
    presolve/Presolve.cpp
    presolve/PresolveComponent.cpp
    presolve/PresolveAnalysis.cpp
    presolve/PresolveUtils.cpp
    presolve/HPreData.cpp
    simplex/HCrash.cpp
    simplex/HDual.cpp
    simplex/HDualRHS.cpp
    simplex/HDualRow.cpp
    simplex/HDualMulti.cpp
    simplex/HFactor.cpp
    simplex/HFactorDebug.cpp
    simplex/HighsSimplexAnalysis.cpp
    simplex/HighsSimplexInterface.cpp
    simplex/HMatrix.cpp
    simplex/HPrimal.cpp
    simplex/HQPrimal.cpp
    simplex/HSimplex.cpp
    simplex/HSimplexDebug.cpp
    simplex/HSimplexReport.cpp
    simplex/HVector.cpp
    test/DevKkt.cpp
    test/KktCh2.cpp
    util/HighsSort.cpp
    util/HighsUtils.cpp
    util/HighsMatrixPic.cpp
    util/stringutil.cpp
    interfaces/highs_c_api.cpp)

set(headers
    ../external/filereaderlp/builder.hpp
    ../external/filereaderlp/model.hpp
    ../external/filereaderlp/reader.hpp
    ipm/IpxStatus.h
    io/Filereader.h
    io/FilereaderLp.h
    io/FilereaderEms.h
    io/FilereaderMps.h
    io/HMpsFF.h
    io/HMPSIO.h
    io/HighsIO.h
    io/LoadOptions.h
    lp_data/HConst.h
    lp_data/HStruct.h
    lp_data/HighsAnalysis.h
    lp_data/HighsDebug.h
    lp_data/HighsInfo.h
    lp_data/HighsLp.h
    lp_data/HighsLpUtils.h
    lp_data/HighsModelUtils.h
    lp_data/HighsModelObject.h
    lp_data/HighsModelObjectUtils.h
    lp_data/HighsOptions.h
    lp_data/HighsRanging.h
    lp_data/HighsSolution.h
    lp_data/HighsSolutionDebug.h
    lp_data/HighsSolve.h
    lp_data/HighsStatus.h
    mip/HighsMipSolver.h
    mip/SolveMip.h
    simplex/HApp.h
    simplex/FactorTimer.h
    simplex/HCrash.h
    simplex/HDual.h
    simplex/HDualRow.h
    simplex/HDualRHS.h
    simplex/HFactor.h
    simplex/HFactorDebug.h
    simplex/HighsSimplexAnalysis.h
    simplex/HighsSimplexInterface.h
    simplex/HMatrix.h
    simplex/HPrimal.h
    simplex/HQPrimal.h
    simplex/HSimplex.h
    simplex/HSimplexReport.h
    simplex/HSimplexDebug.h
    simplex/SimplexConst.h
    simplex/SimplexStruct.h
    presolve/Presolve.h
    presolve/PresolveComponent.h
    presolve/PresolveAnalysis.h
    presolve/PresolveUtils.h
    presolve/HPreData.h
    simplex/HVector.h
    test/DevKkt.h
    test/KktCh2.h
    util/stringutil.h
    util/HighsComponent.h
    util/HighsRandom.h
    util/HighsSort.h
    util/HighsTimer.h
    util/HighsUtils.h
    util/HighsMatrixPic.h
    Highs.h
    interfaces/highs_c_api.h
)

if (GAMS_FOUND)
    set(sources ${sources} interfaces/GAMS.cpp)
    set(sources ${sources} "${GAMS_ROOT}/apifiles/C/api/gmomcc.c")
    set(sources ${sources} "${GAMS_ROOT}/apifiles/C/api/gevmcc.c")
    include_directories(${GAMS_ROOT}/apifiles/C/api)
endif()

if (IPX_ON)
    set(headers ${headers} ipm/IpxWrapper.h ${basiclu_headers}
        ${ipx_headers})
    add_library(libbasiclu ${basiclu_sources})
    add_library(libipx ${ipx_sources})
    target_link_libraries(libipx libbasiclu)
else ()
    set(headers ${headers} ipm/IpxWrapperEmpty.h)
endif()

add_library(libhighs ${sources})

if (IPX_ON)
    target_link_libraries(libhighs libipx)
endif()

if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(libhighs PROPERTIES
        VERSION 
        ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
        SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
else()
    # create static highs library with pic
    set_target_properties(libhighs PROPERTIES
        POSITION_INDEPENDENT_CODE on)
endif()

# on UNIX system the 'lib' prefix is automatically added
set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if (IPX_ON)
    set_target_properties(libipx PROPERTIES
        OUTPUT_NAME "ipx"
        MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    set_target_properties(libbasiclu PROPERTIES
        OUTPUT_NAME "basiclu"
        MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

include(CheckLanguage)
check_language("Fortran")
if(FORTRAN_FOUND)
    set(fortransources interfaces/highs_lp_solver.f90)
    set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
    add_library(FortranHighs interfaces/highs_lp_solver.f90)
    target_link_libraries(FortranHighs PUBLIC libhighs)
    install(TARGETS FortranHighs
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        MODULES DESTINATION modules)
    set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/lib")
endif(FORTRAN_FOUND)

check_language("CSharp")
if(CSHARP_FOUND)
    message(STATUS "CSharp supported")
    set(csharpsources
    interfaces/highs_csharp_api.cs)
    add_library(HighsCsharp interfaces/highs_csharp_api.cs)
    target_compile_options(HighsCsharp PUBLIC "/unsafe")
    add_executable(csharpexample ../examples/call_highs_from_csharp.cs)
    target_compile_options(csharpexample PUBLIC "/unsafe")
    target_link_libraries(csharpexample PUBLIC HighsCsharp)
else()
    message(STATUS "No CSharp support")
endif()

if (OSI_FOUND)
    add_library(OsiHighs interfaces/OsiHiGHSSolverInterface.cpp)
    set(headers ${headers} interfaces/OsiHiGHSSolverInterface.hpp)

    target_include_directories(OsiHighs PUBLIC ${OSI_INCLUDE_DIRS})
    target_link_libraries(OsiHighs PUBLIC libhighs ${OSI_LIBRARIES})
    target_compile_options(OsiHighs PUBLIC ${OSI_CFLAGS_OTHER})

    if(${BUILD_SHARED_LIBS})
        set_target_properties(OsiHighs PROPERTIES
           VERSION 
           ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
           SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
    else()
        set_target_properties(OsiHighs PROPERTIES POSITION_INDEPENDENT_CODE on)
    endif()

    install(TARGETS OsiHighs
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

    set_target_properties(OsiHighs PROPERTIES INSTALL_RPATH
        "${CMAKE_INSTALL_PREFIX}/lib")

    #configure and install the pkg-config file
    configure_file(${HIGHS_SOURCE_DIR}/osi-highs.pc.in
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/osi-highs.pc" @ONLY)
    install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/osi-highs.pc"
        DESTINATION lib/pkg-config)
endif()


# set the install rpath to the installed destination
set_target_properties(highs PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/lib")

# install the header files of highs
foreach ( file ${headers} )
    get_filename_component( dir ${file} DIRECTORY )
    install( FILES ${file} DESTINATION include/${dir} )
endforeach()
install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION include)

if (IPX_ON)
    install(TARGETS libhighs libipx libbasiclu EXPORT highs-targets
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            INCLUDES DESTINATION include)

    # Add library targets to the build-tree export set
    export(TARGETS libhighs libipx libbasiclu
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
else()
    install(TARGETS libhighs EXPORT highs-targets
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            INCLUDES DESTINATION include)

    # Add library targets to the build-tree export set
    export(TARGETS libhighs 
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
endif()


# Configure the config file for the build tree:
# Either list all the src/* directories here, or put explicit paths in all the
# include statements.
# M reckons that the latter is more transparent, and I'm inclined to agree.
set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

# Configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

# Configure the pkg-config file for the install
configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

# Install the targets of the highs export group, the config file so that other
# cmake-projects can link easily against highs, and the pkg-config flie so that
# other projects can easily build against highs
install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION 
    lib/cmake/highs)
install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" 
    DESTINATION lib/cmake/highs)
install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" 
    DESTINATION lib/pkg-config)

else() 

# FAST_BUILD is set to on.
# At the moment used only for gradually updating the CMake targets build and install / export.
# Define library in modern CMake using target_*()
# No interfaces (apart from c); No ipx; New (short) ctest instances.
add_library(libhighs)
target_sources(libhighs PRIVATE 
    ../external/filereaderlp/reader.cpp
    io/Filereader.cpp
    io/FilereaderLp.cpp
    io/FilereaderEms.cpp
    io/FilereaderMps.cpp
    io/HighsIO.cpp
    io/HMPSIO.cpp
    io/HMpsFF.cpp
    io/LoadOptions.cpp
    lp_data/Highs.cpp
    lp_data/HighsDebug.cpp
    lp_data/HighsInfo.cpp
    lp_data/HighsLp.cpp
    lp_data/HighsLpUtils.cpp
    lp_data/HighsModelUtils.cpp
    lp_data/HighsRanging.cpp
    lp_data/HighsSolution.cpp
    lp_data/HighsSolutionDebug.cpp
    lp_data/HighsSolve.cpp
    lp_data/HighsStatus.cpp
    lp_data/HighsOptions.cpp
    mip/HighsMipSolver.cpp
    mip/SolveMip.cpp
    presolve/Presolve.cpp
    presolve/PresolveComponent.cpp
    presolve/PresolveAnalysis.cpp
    presolve/PresolveUtils.cpp
    presolve/HPreData.cpp
    simplex/HCrash.cpp
    simplex/HDual.cpp
    simplex/HDualRHS.cpp
    simplex/HDualRow.cpp
    simplex/HDualMulti.cpp
    simplex/HFactor.cpp
    simplex/HFactorDebug.cpp
    simplex/HighsSimplexAnalysis.cpp
    simplex/HighsSimplexInterface.cpp
    simplex/HMatrix.cpp
    simplex/HPrimal.cpp
    simplex/HQPrimal.cpp
    simplex/HSimplex.cpp
    simplex/HSimplexDebug.cpp
    simplex/HSimplexReport.cpp
    simplex/HVector.cpp
    test/KktCh2.cpp
    test/DevKkt.cpp
    util/HighsSort.cpp
    util/HighsUtils.cpp
    util/HighsMatrixPic.cpp
    util/stringutil.cpp
    interfaces/highs_c_api.cpp)

target_include_directories(libhighs PUBLIC . ipm/ io/ lp_data/ mip/ simplex/ 
                           presolve/ util/ test/ interfaces/ ../external/)
target_include_directories(libhighs PUBLIC ${HIGHS_BINARY_DIR})

if (${HiGHSRELEASE} ON)
    target_compile_options(libhighs PRIVATE "-O3" "-DNDEBUG")
endif()

# target_compile_options(libhighs PRIVATE "-Wall")
# target_compile_options(libhighs PRIVATE "-Wunused")

set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

endif()
