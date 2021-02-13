# The top-level CMake file is there to bring all modules into scope. That
# means, adding the subdirectories for all CMake projects in this tree, and
# finding external libraries and turning them into imported targets.

cmake_minimum_required(VERSION 3.15)

project(HIGHS VERSION 1.0 LANGUAGES CXX C)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). 
    Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Use OpenMP if version is high enough.
option(OPENMP "OpenMP" OFF)
find_package(OpenMP 3.0)
if(OpenMP_CXX_FOUND)
    message(STATUS "Linking with OpenMP... ")
    set(OPENMP ON)
endif()

# Fast build: No interfaces (apart from c); No ipx; New (short) ctest instances, 
# static library and exe without PIC. Used for gradually updating the CMake 
# targets build and install / export.

# Note: Not ON by default. Still need to fix the exe generation dir: rename from 
# app/ to bin/
option(EXP "Experimental mode: run unit tests with doctest." OFF)
option(FAST_BUILD "Fast build: only build static lib and exe and quick test." OFF)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/lib)

# set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

if (NOT FAST_BUILD)

set(HIGHS_VERSION_PATCH 0)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

find_package(PkgConfig)

include(CPack)
set(CPACK_RESOURCE_FILE_LICENSE "${HIGHS_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/lib)

option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "Export all symbols into the DLL" ON)

find_program(GIT git)
if((GIT) AND (EXISTS ${HIGHS_SOURCE_DIR}/.git))
    execute_process(
        COMMAND ${GIT} describe --always --dirty
        WORKING_DIRECTORY ${HIGHS_SOURCE_DIR}
        OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^.*-g" "" GITHASH ${GITHASH})
else()
    set(GITHASH "n/a")
endif()
message(STATUS "Git hash: " ${GITHASH})

string(TIMESTAMP TODAY "%Y-%m-%d")
message(STATUS "Compilation date: " ${TODAY})

set(EXT_PRESOLVE_DIR "" CACHE STRING "Path to the external presolving library")
if(EXT_PRESOLVE_DIR)
    message(STATUS "External presolve: " ${EXT_PRESOLVE_DIR})
    include_directories(${EXT_PRESOLVE_DIR}/src)
    include_directories(${EXT_PRESOLVE_DIR}/external)
endif()



# For the moment keep above coverage part in case we are testing at CI.
option(CI "CI extended tests" ON)

# Coverage part
# 'make coverage' to start the coverage process
option(HIGHS_COVERAGE "Activate the code coverage compilation" OFF)
if (HIGHS_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}    -O0 --coverage")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG}  -O0 --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O0 --coverage")
  endif ()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    message(FATAL_ERROR "Warning: to enable coverage, you must compile in DEBUG mode")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (WIN32)
    message(FATAL_ERROR "Error: code coverage analysis is only available under Linux for now.")
  endif ()
  
  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)

  if (NOT GCOV_PATH)
    message(FATAL_ERROR "gcov not found! Please install lcov and gcov. Aborting...")
  endif ()
  
  if (NOT LCOV_PATH)
    message(FATAL_ERROR "lcov not found! Please install lcov and gcov. Aborting...")
  endif ()
  
  if (NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found! Please install lcov and gcov. Aborting...")
  endif ()
  
  # Capturing lcov counters and generating report
  if (NOT CI)
    add_custom_target(coverage
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/src --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${LCOV_PATH} --remove "*/usr/include/*" --output-file ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${GENHTML_PATH} -o ${CMAKE_BINARY_DIR}/coverage ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/coverage.info.cleaned
            VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.
    Processing code coverage counters and generating report.
    You can zip the directory ${CMAKE_BINARY_DIR}/coverage and upload the content to a web server.")
  else()
    add_custom_target(ci_cov
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/src --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
		    VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.")
  endif()
endif ()

option(SCIP_DEV "SCIP  development" OFF)
if (NOT (${SCIP_DEV} STREQUAL  "OFF"))
    message(STATUS "SCIP  development: " ${SCIP_DEV})
endif()

option(HiGHSDEV "HiGHS development" OFF)
if (NOT (${HiGHSDEV} STREQUAL  "OFF"))
    message(STATUS "HiGHS development: " ${HiGHSDEV})
endif()

set(OSI_ROOT "" CACHE PATH "Osi root folder.")
if (NOT "${OSI_ROOT}" STREQUAL "")
    # if OSI_ROOT is set, then overwrite PKG_CONFIG_PATH
    message(STATUS "OSI root folder set: " ${OSI_ROOT})
    set(ENV{PKG_CONFIG_PATH}  "${OSI_ROOT}/lib/pkgconfig")
endif ()
unset(OSI_ROOT CACHE)

pkg_check_modules(OSI osi)
if (OSI_FOUND)
    # need to come before adding any targets (add_executable, add_library)
    link_directories(${OSI_LIBRARY_DIRS})
    include_directories(${OSITEST_INCLUDE_DIRS})
endif (OSI_FOUND)

set(IPX_ON ON)
#option(IPX_ON OFF)
#if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux" AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5))
#    set(IPX_ON ON)
#endif()

include(CheckLanguage)
check_language("Fortran")
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    set(FORTRAN_FOUND ON)
else()
    set(FORTRAN_FOUND OFF)
endif(CMAKE_Fortran_COMPILER)

check_language("CSharp")
if(CMAKE_CSharp_COMPILER)
    enable_language(CSharp)
    set(CSHARP_FOUND ON)
else()
    set(CSHARP_FOUND OFF)
endif(CMAKE_CSharp_COMPILER)

option(GAMS_ROOT "GAMS root folder." "OFF")
if (NOT (${GAMS_ROOT} STREQUAL "OFF"))
    message(STATUS "GAMS root folder set: " ${GAMS_ROOT})
    set(GAMS_FOUND ON)
else ()
    set(GAMS_FOUND OFF)
endif ()

# whether to use shared or static libraries
option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(HiGHSRELEASE ON)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)
include_directories(
    ${HIGHS_BINARY_DIR}
    ${HIGHS_SOURCE_DIR}/app
    ${HIGHS_SOURCE_DIR}/external
    ${HIGHS_SOURCE_DIR}/src
    ${HIGHS_SOURCE_DIR}/src/io
    ${HIGHS_SOURCE_DIR}/src/ipm/ipx
    ${HIGHS_SOURCE_DIR}/src/ipm/basiclu
    ${HIGHS_SOURCE_DIR}/src/lp_data
    ${HIGHS_SOURCE_DIR}/src/mip
    ${HIGHS_SOURCE_DIR}/src/presolve
    ${HIGHS_SOURCE_DIR}/src/simplex
    ${HIGHS_SOURCE_DIR}/src/test
    ${HIGHS_SOURCE_DIR}/src/util)


# explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()

# If Visual Studio targets are being built.
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Function to set compiler flags on and off easily.
include(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

# usage: turn pedantic on for even more warnings.
enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-Wno-unused-parameter")
enable_cxx_compiler_flag_if_supported("-pedantic")

#if(CMAKE_BUILD_TYPE STREQUAL Debug OR CMAKE_BUILD_TYPE STREQUAL debug)
#    enable_cxx_compiler_flag_if_supported("-D_GLIBCXX_DEBUG")
#endif()

if (IPX_ON)
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Targets
enable_testing()
add_subdirectory(app)
add_subdirectory(check)
add_subdirectory(src)


else()
message(STATUS "FAST BUILD set to on.
 Note: The HiGHS team is preparing for our first official release. If you
       experience any issues please let us know via email or on GitHub.")

# todo: remove patch once release is official. 
# add here to avoid compile errors.
set(HIGHS_VERSION_PATCH 0)
configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)

# condition added to src/CMakeLists.txt
add_subdirectory(src)

# No changes in app/ apart from a relative path
add_subdirectory(app)

# check/ not added here, instead define fewer tests: 
# build, 3 feas, 1 infeas, 1 unbounded. 1 parallel, 1 no presolve
enable_testing()

# Check whether targets build OK.
add_test(NAME highs-lib-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target libhighs
         )

set_tests_properties(highs-lib-build
                     PROPERTIES
                     RESOURCE_LOCK libhighs)

add_test(NAME highs-exe-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target highs
         )

set_tests_properties(highs-exe-build
                     PROPERTIES
                     RESOURCE_LOCK highs)

set(successInstances
    "25fv47\;2888\; 5.5018458883\;"
    "80bau3b\;3760\; 9.8722419241\;"
    "greenbea\;5249\;-7.2555248130\;")

set(optionsInstances
    "adlittle\;81\; 2.2549496316\;")

set(infeasibleInstances
    "bgetam\;        infeasible")

set(unboundedInstances
     "gas11\;         unbounded")
     
# define settings
set(settings
    ""
    "--presolve=off"
    "--parallel=on")

# define function to add tests
# More Modern CMake: avoid macros if you can
function(add_instance_tests instances solutionstatus setting)
# loop over the instances
foreach(instance ${${instances}})
    # add default tests
    # treat the instance as a tuple (list) of two values
    list(GET instance 0 name)
    list(GET instance 1 iter)

    if(${solutionstatus} STREQUAL "Optimal")
        list(GET instance 2 optval)
    endif()

    # specify the instance and the settings load command
    set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps")

    add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
            ${inst})

    set_tests_properties (${name}${setting} PROPERTIES
            DEPENDS unit_tests_all)
    set_tests_properties (${name}${setting} PROPERTIES
            PASS_REGULAR_EXPRESSION
            "Model   status      : ${solutionstatus}")

    if(${solutionstatus} STREQUAL "Optimal")
        if(${setting} STREQUAL "--presolve=off")
            set_tests_properties (${name}${setting} PROPERTIES
                    PASS_REGULAR_EXPRESSION
                    "Simplex   iterations: ${iter}\nObjective value     : ${optval}")
        else()
            set_tests_properties (${name}${setting} PROPERTIES
                    PASS_REGULAR_EXPRESSION
                    "Objective value     : ${optval}")
        endif()
    endif()
endforeach(instance)

endfunction()

# add tests for success and fail instances
add_instance_tests(successInstances "Optimal" " ")
add_instance_tests(failInstances "Fail" " ")
add_instance_tests(infeasibleInstances "Infeasible" " ")
add_instance_tests(unboundedInstances "Unbounded" " ")

foreach(setting ${settings})
    add_instance_tests(optionsInstances "Optimal" ${setting})
endforeach()

if (EXP) 
    add_executable(doctest)
    #target_sources(doctest PRIVATE check/doctest/TestPresolveColumnSingletons.cpp)
    target_sources(doctest PRIVATE check/doctest/TestPresolveIssue.cpp)
    target_include_directories(doctest PRIVATE external)
    target_link_libraries(doctest libhighs)
endif()

# Dashboards later on.
# include(CTest)

endif()

# Modern CMake link in FAST_BUILD mode
# All uses of target_link_libraries with a target must be either 
# all-keyword or all-plain.
if (OPENMP)
    if (FAST_BUILD)
        target_link_libraries(libhighs PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_link_libraries(libhighs OpenMP::OpenMP_CXX)
    endif()
endif()

# # Comment out for scaffold/ tests
# add_subdirectory(scaffold)

# # Only needed if it defines a target like a main executable. For methods stick
# # to header-only.
# add_subdirectory(dev_presolve)

