# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)
option(BUILD_SHARED_EXTRAS_LIB "Build shared extras library" ON)

# configure the version
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
include(vcpkg)
include(sources)
include(set-version)

set_version(VERSION ${CMAKE_CURRENT_LIST_DIR}/..)
project(highs_extras VERSION ${VERSION} LANGUAGES C CXX)

set(HIGHS_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(HIGHS_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(HIGHS_VERSION_PATCH ${PROJECT_VERSION_PATCH})
configure_file(${CMAKE_CURRENT_LIST_DIR}/../highs/HConfig.h.in ${PROJECT_BINARY_DIR}/HConfig.h)
message(STATUS "Configuring highs_extras version ${VERSION}")

set(highs_extras_sources
    HighsExtrasApi.cpp
    HighsExtrasExternalDeps.cpp)

set(highs_extras_headers
    HighsExtrasApi.h
    HighsExtrasApiBinding.h
    HighsExtrasExternalDeps.h
    HighsMetis.h
    OrderingPrint.h)

# add optional dependencies (source code)
if (HIPO)
  include(FindHipoDeps)
  highs_configure_blas()

  set(highs_extras_sources ${highs_extras_sources} 
    ${amd_sources} ${blas_sources} ${metis_sources} ${rcm_sources})

  set(highs_extras_headers ${highs_extras_headers} 
    ${amd_headers} ${blas_headers} ${metis_headers} ${rcm_headers})
endif()
# end optional dependencies (source code)

if (BUILD_SHARED_EXTRAS_LIB)
  add_library(highs_extras SHARED
      ${highs_extras_sources}
      ${highs_extras_headers})

  # build flag to indicate shared library is being built, so API symbols are exported
  target_compile_definitions(highs_extras PRIVATE HIGHS_EXTRAS_LIBRARY_BUILD)
else()
  add_library(highs_extras STATIC
      ${highs_extras_sources}
      ${highs_extras_headers})
endif()

# use C++11 standard
target_compile_features(highs_extras PUBLIC cxx_std_11)

set_target_properties(highs_extras PROPERTIES
    CXX_EXTENSIONS OFF
    OUTPUT_NAME "highs_extras"
    POSITION_INDEPENDENT_CODE ON)

target_include_directories(highs_extras PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

target_compile_definitions(highs_extras PRIVATE HIGHS_EXTRAS_VERSION="${VERSION}")

if(MSVC)
  target_compile_options(highs_extras PRIVATE "/bigobj")
else()
  target_compile_options(highs_extras PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-ftemplate-depth=2048>)
endif()

# add optional dependencies (linker)
if(HIPO)
  highs_link_blas(highs_extras)
endif()
# end optional dependencies (linker)

if(BUILD_SHARED_EXTRAS_LIB)
  set_target_properties(highs_extras PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)

  if(NOT WIN32)
    if(APPLE)
      set(highs_extras_rpath "@loader_path")
    else()
      set(highs_extras_rpath "$ORIGIN")
    endif()

    set_target_properties(highs_extras PROPERTIES
      BUILD_RPATH "${highs_extras_rpath}"
      INSTALL_RPATH "${highs_extras_rpath}"
      INSTALL_RPATH_USE_LINK_PATH TRUE)
  else()
    # no "lib" prefix with MinGW: the loader looks for highs_extras.dll
    set_target_properties(highs_extras PROPERTIES
      PREFIX ""
      WINDOWS_EXPORT_ALL_SYMBOLS OFF)
  endif()
endif()

# This CMake file is only inlcluded if not (HIPO_PYTHON)
install(TARGETS highs_extras
    EXPORT highs-targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT libs)

if(DEFINED HIGHS_BINARY_DIR)
  export(TARGETS highs_extras
    APPEND FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
  if(HIPO)
    highs_export_blas_dependency()
    highs_export_blas_pkgconfig()
  endif()
endif()