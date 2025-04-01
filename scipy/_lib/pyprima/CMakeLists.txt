cmake_minimum_required (VERSION 3.13)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (prima Fortran)

option (BUILD_SHARED_LIBS "shared/static" ON)

include (GNUInstallDirs)

# the compiler enables an executable stack because of nested functions and this is fine
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  include (CheckLinkerFlag)
  check_linker_flag (Fortran "-Wl,--no-warn-execstack" HAVE_WARN_EXECSTACK)
endif ()

# Set additional Fortran compiler flags
# 0. See https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html for compiler IDs.
# 1. We require the compilers to allocate arrays on the heap instead of the stack, which is
# slower (does not matter for DFO applications) but can avoid memory errors on large problems.
# 2. We require the compilers to compile the solvers so that they can be called recursively.
# See https://fortran-lang.discourse.group/t/frecursive-assume-recursion-and-recursion-thread-safety
option (PRIMA_HEAP_ARRAYS "allocate arrays on heap" ON)
if (PRIMA_HEAP_ARRAYS)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")  # gfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-stack-arrays -frecursive")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")  # Intel compilers
    if (WIN32)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /heap-arrays /assume:recursion")
    else ()
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays -assume recursion")
    endif ()
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")  # nagfor
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive")  # What about stack/heap?
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")  # flang-new
    # See https://github.com/llvm/llvm-project/issues/88344
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-stack-arrays -mmlir -fdynamic-heap-array")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")  # Classic Flang and AOCC Flang
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "ARMClang")  # ARM Flang
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-stack-arrays -Mrecursive")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")  # nvfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mnostack_arrays -Mrecursive")
  endif ()
endif ()

# Set additional linker flags
# Zaikun 20240217: Fix https://github.com/libprima/prima/issues/158, which is due to the new linker
# implemented in Xcode 15 on macOS. It happens only if the Fortran compiler is ifort.
# An alternative is `add_link_options("-ld_classic")`, which forces Xcode to use the old linker.
# Will CMake adapt itself to the new linker later? Will a newer version of CMake make the fix unnecessary?
# See
# https://developer.apple.com/documentation/xcode-release-notes/xcode-15-release-notes#Linking
# https://stackoverflow.com/questions/77525544/apple-linker-warning-ld-warning-undefined-error-is-deprecated
# https://medium.com/@hackingwithcode/cmake-and-xcode-15-solving-the-undefined-error-puzzle-3c847e6d1008
# Zaikun 20240501: "undefined,dynamic_lookup" seems included by default since CMake 3.28.1. We keep
# the following code in case it is not included for some versions.
if((APPLE) AND (CMAKE_Fortran_COMPILER_ID MATCHES "Intel"))
    add_link_options("-Wl,-undefined,dynamic_lookup")
endif()
# Zaikun 20240501: Fix "ld: Assertion failed: (resultIndex < sectData.atoms.size())", which happens
# when building Python wheels on macOS 13 and 14. This is a bug of Xcode 15.0 on macOS and fixed in
# Xcode 15.1 beta according to https://github.com/Homebrew/homebrew-core/issues/145991.
# See also
# https://forums.developer.apple.com/forums/thread/737707
# https://github.com/hansec/OpenFUSIONToolkit/pull/29
if(APPLE)
    # Get the Xcode version
    execute_process(
        COMMAND xcodebuild -version
        OUTPUT_VARIABLE XCODE_OUTPUT
    )
    string(REGEX MATCH "Xcode ([0-9]+\\.[0-9]+(\\.[0-9]+)?)" XCODE_VERSION "${XCODE_OUTPUT}")
    # Set linker flags for Xcode 15.0 or 15.0.x
    if(XCODE_VERSION MATCHES "Xcode 15\\.0(\\.[0-9]+)?")
        add_link_options("-Wl,-ld_classic")
    endif()
endif()

# For running tests with gdb. $_exitcode != 0 means the program ran without exiting
# normally, and in this case we want to show a stack trace
file(WRITE ${CMAKE_BINARY_DIR}/cmdfile.gdb "init-if-undefined $_exitcode = 0
run
set language c
if $_exitcode != 0
  where
end
quit $_exitcode
")

option(PRIMA_ENABLE_EXAMPLES "build examples by default" OFF)
add_custom_target (examples)
enable_testing ()

option(PRIMA_ENABLE_TESTING "build tests" OFF)
add_custom_target (tests)
add_dependencies(tests examples)
add_subdirectory(fortran)

option (PRIMA_ENABLE_C "C binding" ON)

if (PRIMA_ENABLE_C)
  enable_language(C)
  add_subdirectory(c)
  set(primac_target "primac")
endif ()

# Get the version number
find_package(Git)
set(IS_REPO FALSE)
if(GIT_EXECUTABLE)
  # --always means git describe will output the commit hash if no tags are found
  # This is usually the case for forked repos since they do not clone tags by default.
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE PRIMA_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE GIT_RESULT ERROR_QUIET)
  if (GIT_RESULT EQUAL 0)
    set(IS_REPO TRUE)
  endif()
endif()
if(NOT GIT_EXECUTABLE OR NOT IS_REPO)
  # If git is not available, or this isn't a repo, that may indicate we are building
  # on macports which downloads the bundle from github (which uses git archive) and
  # so the version number should be in .git-archival.txt.
  # Alternatively it might mean that we're building the Python bindings, in which case
  # the version is output in _version.txt. I know, it's complicated. I don't make the rules.
  if(EXISTS _version.txt)
    file(STRINGS _version.txt PRIMA_VERSION)
  else()
    file(STRINGS .git-archival.txt PRIMA_VERSION)
    if(PRIMA_VERSION MATCHES "describe")
      message(WARNING "No git detected and .git-archival.txt does not contain a version number")
      set(PRIMA_VERSION "unknown")
    endif()
  endif()

endif()
# Remove the leading v from PRIMA_VERSION, if it contains one.
string(REGEX REPLACE "^v" "" PRIMA_VERSION ${PRIMA_VERSION})
message(STATUS "Setting PRIMA version to ${PRIMA_VERSION}")

option (PRIMA_ENABLE_PYTHON "Python binding" OFF)
if (PRIMA_ENABLE_PYTHON)
  if(NOT PRIMA_ENABLE_C)
    message(FATAL_ERROR "Building Python bindings requires C bindings. Please turn on PRIMA_ENABLE_C")
  endif()
  if(BUILD_SHARED_LIBS)
    # This will include libprimaf, libprimafc, and libprimac into the compiled Python binding, removing the need
    # to properly set the rpath or find those libraries at runtime.
    # Even if we did make it successfully build with shared libraries, delocate/auditwheel will copy them into the
    # bindings anyway
    message(FATAL_ERROR "Building Python bindings requires static libraries. Please disable BUILD_SHARED_LIBS")
  endif()
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(WARNING "Compiling Python bindings with compilers other than GNU has not been tested and no support is planned at this time")
  endif()
  enable_language(CXX)
  add_subdirectory(python)
endif ()

install(
    TARGETS primaf ${primac_target}
    EXPORT prima-targets
    INCLUDES DESTINATION include
)

install(
    EXPORT prima-targets
    FILE prima-targets.cmake
    NAMESPACE prima::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/prima-config.cmake.in
    ${CMAKE_BINARY_DIR}/prima-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/prima-config-version.cmake
    VERSION ${PRIMA_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/prima-config.cmake
        ${CMAKE_BINARY_DIR}/prima-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/prima
)
