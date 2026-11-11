if (NOT BUILD_CXX)
  return()
endif()

# Define library.
include(sources)
set(sources ${highs_sources} ${cupdlp_sources} ${ipx_sources} ${basiclu_sources} ${hipo_sources} ${factor_highs_sources} ${hipo_util_sources})
set(headers ${highs_headers} ${cupdlp_headers} ${ipx_headers} ${basiclu_headers} ${hipo_headers} ${factor_highs_headers} ${hipo_util_headers})

# Configure the config windows version file
if(MSVC)
  string(REPLACE "." ","  PROJECT_RC_VERSION "${PROJECT_VERSION}")
  configure_file(${PROJECT_SOURCE_DIR}/version.rc.in
    "${HIGHS_BINARY_DIR}/version.rc" @ONLY)
  set(win_version_file ${HIGHS_BINARY_DIR}/version.rc)
else()
  set(win_version_file)
endif()

# Outdated CMake approach: update in progress
if(NOT FAST_BUILD)

  add_library(libhighs ${sources} ${headers} ${win_version_file})
  target_include_directories(libhighs PRIVATE ${include_dirs})

  if(APPLE)
    # linker was taking ~1 hour, unity build reduces it to ~1 minute
    set_target_properties(libhighs PROPERTIES
      UNITY_BUILD ON
      UNITY_BUILD_BATCH_SIZE 8)

    if(NOT CMAKE_VERSION VERSION_LESS 3.20)
      set_target_properties(libhighs PROPERTIES UNITY_BUILD_UNIQUE_ID "HIGHS_UNITY_ID")
    endif()
  endif()

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(libhighs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
    if(MINGW)
      target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
    endif()
  else()
    # create static highs library with pic
    set_target_properties(libhighs PROPERTIES
      POSITION_INDEPENDENT_CODE on)
    target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
  endif()

  # on UNIX system the 'lib' prefix is automatically added
  set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    PDB_NAME "libhighs"
  )

  if(ZLIB AND ZLIB_FOUND)
    target_link_libraries(libhighs ZLIB::ZLIB)
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(ZLIB)")
  endif()

  # set the install rpath to the installed destination
  if(APPLE)
      set_target_properties(libhighs PROPERTIES INSTALL_RPATH "@loader_path")
  elseif (UNIX)
      set_target_properties(libhighs PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()

  # install the header files of highs
  foreach(file ${headers})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  if(UNIX)
    # target_compile_options(libhighs PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libhighs PRIVATE "-Wno-return-type-c-linkage")
    target_compile_options(libhighs PRIVATE "-Wno-return-type" "-Wno-switch")

    target_compile_options(libhighs PRIVATE "-Wno-unused-variable")
    target_compile_options(libhighs PRIVATE "-Wno-unused-const-variable")

    # target_compile_options(libhighs PRIVATE "-Wno-sign-compare")
    # target_compile_options(libhighs PRIVATE "-Wno-logical-op-parentheses")

    # target_compile_options(libipx PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libipx PRIVATE "-Wno-return-type-c-linkage")
    # target_compile_options(libipx PRIVATE "-Wno-return-type" "-Wno-switch")

    # target_compile_options(libipx PRIVATE "-Wno-sign-compare")
    # target_compile_options(libipx PRIVATE "-Wno-logical-op-parentheses")
  endif()

  install(TARGETS libhighs EXPORT highs-targets
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # Add library targets to the build-tree export set
  export(TARGETS libhighs
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")

  # Configure the config file for the build tree:
  # Either list all the highs/* directories here, or put explicit paths in all the
  # include statements.
  # M reckons that the latter is more transparent, and I'm inclined to agree.
  set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/highs" "${HIGHS_BINARY_DIR}")
  configure_file(${PROJECT_SOURCE_DIR}/cmake/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

  # Configure the config file for the install
  set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${CMAKE_INSTALL_INCLUDEDIR}/highs")
  configure_file(${PROJECT_SOURCE_DIR}/cmake/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

  # Configure the pkg-config file for the install
  configure_file(${PROJECT_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

  # Install the targets of the highs export group, the config file so that other
  # cmake-projects can link easily against highs, and the pkg-config file so that
  # other projects can easily build against highs
  install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

else()
  # FAST_BUILD is set to on.
  # At the moment used only for gradually updating the CMake targets build and
  # install / export.
  # Define library in modern CMake using target_*()
  # No interfaces (apart from c); No ipx; New (short) ctest instances.
  add_library(highs)
  add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

  if(APPLE)
    # linker was taking ~1 hour, unity build reduces it to ~1 minute
    set_target_properties(highs PROPERTIES
      UNITY_BUILD ON
      UNITY_BUILD_BATCH_SIZE 8)

    if(NOT CMAKE_VERSION VERSION_LESS 3.20)
      set_target_properties(highs PROPERTIES UNITY_BUILD_UNIQUE_ID "HIGHS_UNITY_ID")
    endif()
  endif()

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(highs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
  endif()

  set_target_properties(highs PROPERTIES POSITION_INDEPENDENT_CODE ON)

  if(APPLE)
      set_target_properties(highs PROPERTIES INSTALL_RPATH "@loader_path")
  elseif (UNIX)
      set_target_properties(highs PROPERTIES INSTALL_RPATH "$ORIGIN")
  endif()

  if (BUILD_SHARED_EXTRAS_LIB)
    if (NOT HIPO_PYTHON)
      if(APPLE)
        set_target_properties(highs PROPERTIES
          INSTALL_RPATH "@loader_path;@loader_path/../${CMAKE_INSTALL_LIBDIR}"
          BUILD_RPATH   "$<TARGET_FILE_DIR:highs_extras>"
        )
      elseif(UNIX)
        set_target_properties(highs PROPERTIES
          INSTALL_RPATH "$ORIGIN"
          BUILD_RPATH   "$<TARGET_FILE_DIR:highs_extras>"
        )
      endif()
    endif()
  else()
    # PUBLIC: propagate link dependency in the highs library, so that 
    # consumers don't need to explicitly link against highs_extras
    target_link_libraries(highs PUBLIC highs_extras)
  endif()

  target_sources(highs PRIVATE ${sources} ${headers} ${win_version_file})

  # Optional GPU (CUDA or HIP)
  if (CUPDLP_GPU OR HIPDLP_GPU)

    if (HIPDLP_HIP)
      # HIP backend — pdhg.cu is compiled as HIP source.
      # HIPDLP_GPU activates GPU code in HiPDLP headers.
      # CUPDLP_CPU keeps cuPDLP-C in CPU mode so it doesn't pull CUDA headers.
      target_compile_definitions(highs PRIVATE HIPDLP_GPU USE_HIP CUPDLP_CPU)

      set_source_files_properties(pdlp/hipdlp/pdhg.cu PROPERTIES LANGUAGE HIP)

      target_sources(highs PRIVATE pdlp/hipdlp/pdhg.cu)

      set_target_properties(highs PROPERTIES
        HIP_SEPARABLE_COMPILATION ON)

      if (WIN32)
        target_link_libraries(highs PRIVATE ${GPU_LIBRARY})
      else()
        target_link_libraries(highs PRIVATE ${GPU_LIBRARY} m)
      endif()

    else()
      # CUDA backend
      target_include_directories(highs PUBLIC "$<BUILD_INTERFACE:${CMAKE_CUDA_PATH}/include>")

      if (HIGHS_GPU_LIB)
        # like in cuda CMakeLists
        enable_language(CXX CUDA)

        # Remove /MP flag from CUDA targets (nvcc doesn't support it)
        if(MSVC)
          string(REPLACE "/MP" "" CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
          string(REPLACE "/bigobj" "" CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")

          foreach(flag_var
            CMAKE_CUDA_FLAGS_DEBUG
            CMAKE_CUDA_FLAGS_RELEASE
            CMAKE_CUDA_FLAGS_RELWITHDEBINFO
            CMAKE_CUDA_FLAGS_MINSIZEREL)
            string(REPLACE "/MP" "" ${flag_var} "${${flag_var}}")
            string(REPLACE "/bigobj" "" ${flag_var} "${${flag_var}}")
          endforeach()
        endif()

        target_sources(highs PRIVATE
          pdlp/cupdlp/cuda/cupdlp_cuda_kernels.cu
          pdlp/cupdlp/cuda/cupdlp_cuda_kernels.cuh
          pdlp/cupdlp/cuda/cupdlp_cudalinalg.cuh
          pdlp/cupdlp/cuda/cupdlp_cudalinalg.cu
          pdlp/cupdlp/cupdlp_cs.c
          pdlp/hipdlp/pdhg.cu)

        set_target_properties(highs PROPERTIES
          CUDA_SEPARABLE_COMPILATION ON
          #CUDA_ARCHITECTURES "60;70;75;80;86" )
          #Deleted 60;70 due to errors
          #nvcc fatal   : Unsupported gpu architecture 'compute_60'
          #nvcc fatal   : Unsupported gpu architecture 'compute_70'
          CUDA_ARCHITECTURES "75;80;86" )

        target_include_directories(highs PUBLIC "$<BUILD_INTERFACE:${CMAKE_CUDA_PATH}/include>")

        if (WIN32)
          target_link_libraries(highs PRIVATE ${GPU_LIBRARY})
        else()
          target_link_libraries(highs PRIVATE ${GPU_LIBRARY} m)
        endif()

      else()
        set(CUPDLP_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/highs/pdlp/cupdlp/")
        add_subdirectory(pdlp/cupdlp/cuda)

        if (WIN32)
          target_link_libraries(highs PRIVATE cudalin ${GPU_LIBRARY})
        else()
          target_link_libraries(highs PRIVATE cudalin ${GPU_LIBRARY} m)
        endif()

        set_target_properties(highs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
      endif()
    endif()
  endif()

  target_include_directories(highs PRIVATE ${include_dirs})

  target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs/hipo>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs/hipo/auxiliary>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs/hipo/factorhighs>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs/hipo/ipm>
  )

  if(MSVC)
    if (NOT HIGHS_GPU_LIB)
      list(APPEND highs_compile_opts
        "/bigobj" # Allow big object
        "/DNOMINMAX"
        "/DWIN32_LEAN_AND_MEAN=1"
        "/D_CRT_SECURE_NO_WARNINGS"
        "/D_CRT_SECURE_NO_DEPRECATE"
        "/MP" # Build with multiple processes
        "/Zc:preprocessor" # Enable preprocessor conformance mode
        "/fp:precise"
        )
    endif()
    if (CMAKE_BUILD_TYPE STREQUAL Release)
      list(APPEND highs_compile_opts
        "/DNDEBUG")
    endif()
    # MSVC warning suppressions
    if (NOT HIGHS_GPU_LIB)
    list(APPEND highs_compile_opts
      "/wd4005" # 'macro-redefinition'
      "/wd4018" # 'expression' : signed/unsigned mismatch
      "/wd4065" # switch statement contains 'default' but no 'case' labels
      "/wd4068" # 'unknown pragma'
      "/wd4101" # 'identifier' : unreferenced local variable
      "/wd4146" # unary minus operator applied to unsigned type, result still unsigned
      "/wd4200" # nonstandard extension used : zero-sized array in struct/union
      "/wd4244" # 'conversion' conversion from 'type1' to 'type2', possible loss of data
      "/wd4251" # 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'
      "/wd4267" # 'var' : conversion from 'size_t' to 'type', possible loss of data
      "/wd4305" # 'identifier' : truncation from 'type1' to 'type2'
      "/wd4307" # 'operator' : integral constant overflow
      "/wd4309" # 'conversion' : truncation of constant value
      "/wd4334" # 'operator' : result of 32-bit shift implicitly converted to 64 bits (was 64-bit shift intended?)
      "/wd4355" # 'this' : used in base member initializer list
      "/wd4477" # 'fwprintf' : format string '%s' requires an argument of type 'wchar_t *'
      "/wd4506" # no definition for inline function 'function'
      "/wd4715" # function' : not all control paths return a value
      "/wd4800" # 'type' : forcing value to bool 'true' or 'false' (performance warning)
      "/wd4996" # The compiler encountered a deprecated declaration.
      )
    endif()
    target_compile_options(highs PRIVATE ${highs_compile_opts})
  # else()
  #   list(APPEND highs_compile_opts "-fwrapv")
  endif()


  if(ZLIB AND ZLIB_FOUND)
    target_include_directories(highs PRIVATE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/zstr>
    )
    target_link_libraries(highs PRIVATE ZLIB::ZLIB)
  endif()

  # allow consumers of highs to be aware of potential external dependencies
  # without forcing them to explicitly call target_include_directories
  target_include_directories(highs PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern>
  )

  # install the header files of highs
  foreach(file ${headers})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # target_compile_options(highs PRIVATE "-Wall")
  # target_compile_options(highs PRIVATE "-Wunused")

  if (UNIX)
    if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      target_compile_options(highs PRIVATE "-Wall")
      target_compile_options(highs PRIVATE "-Wreturn-type")
      target_compile_options(highs PRIVATE "-Wmissing-declarations")
      target_compile_options(highs PRIVATE "-Wno-unused-variable")
      target_compile_options(highs PRIVATE "-Wno-unused-const-variable")
      target_compile_options(highs PRIVATE "-Wno-unused-but-set-variable")
      target_compile_options(highs PRIVATE "-Wno-comment")
      target_compile_options(highs PRIVATE "-Wno-unused-label")

      if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(highs PRIVATE "-Wno-unused-lambda-capture")
      endif()

    else()
      target_compile_options(highs PRIVATE "-Wno-unused-variable")
      target_compile_options(highs PRIVATE "-Wno-unused-const-variable")
    endif()
  endif()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(highs PUBLIC "-Wno-invalid-offsetof")
  endif()

  check_cxx_compiler_flag("-Wno-invalid-offsetof" COMPILER_SUPPORTS_INVALID_OFFSET)
  if(COMPILER_SUPPORTS_INVALID_OFFSET)
      target_compile_options(highs PUBLIC "$<$<COMPILE_LANGUAGE:CXX>:-Wno-invalid-offsetof>")
  endif()

  if (BUILD_DOTNET)

      # see: https://docs.microsoft.com/en-us/dotnet/core/rid-catalog
      # With MSVC, CMAKE_SYSTEM_PROCESSOR reports the host rather than the -A
      # target platform, so prefer the compiler's target architecture.
      if(MSVC AND CMAKE_CXX_COMPILER_ARCHITECTURE_ID)
        string(TOLOWER "${CMAKE_CXX_COMPILER_ARCHITECTURE_ID}" DOTNET_TARGET_ARCH)
      else()
        string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" DOTNET_TARGET_ARCH)
      endif()

      if(DOTNET_TARGET_ARCH MATCHES "^(aarch64|arm64)")
        set(DOTNET_PLATFORM arm64)
      elseif(WIN32 AND (CMAKE_GENERATOR_PLATFORM STREQUAL "Win32" OR CMAKE_SIZEOF_VOID_P EQUAL 4))
        set(DOTNET_PLATFORM x86)
      else()
        set(DOTNET_PLATFORM x64)
      endif()

      if(APPLE)
        set(DOTNET_RID osx-${DOTNET_PLATFORM})
      elseif(UNIX)
        set(DOTNET_RID linux-${DOTNET_PLATFORM})
      elseif(WIN32)
        set(DOTNET_RID win-${DOTNET_PLATFORM})
      else()
        message(FATAL_ERROR "Unsupported system !")
      endif()
      message(STATUS ".Net RID: ${DOTNET_RID}")

      set(DOTNET_PROJECT Highs.Native)
      set(DOTNET_PROJECT_DIR ${PROJECT_BINARY_DIR}/dotnet/${DOTNET_PROJECT})
      file(MAKE_DIRECTORY ${DOTNET_PROJECT_DIR}/runtimes/${DOTNET_RID}/native)

      if (APPLE)
        set(TARGET_FILE_NAME "highs.dylib")
      elseif(UNIX)
        set(TARGET_FILE_NAME "highs.so")
      elseif(WIN32)
        set(TARGET_FILE_NAME "highs.dll")
      endif()

      add_custom_command(TARGET highs POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
          "$<TARGET_FILE:highs>"
        ${DOTNET_PROJECT_DIR}/runtimes/${DOTNET_RID}/native/${TARGET_FILE_NAME}
        COMMENT "Copying to output directory")
  endif()

endif()

if(FORTRAN_FOUND)
  set(fortransources interfaces/highs_fortran_api.f90)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
  add_library(FortranHighs interfaces/highs_fortran_api.f90)

  if(NOT FAST_BUILD)
    target_link_libraries(FortranHighs PUBLIC libhighs)
  else()
    target_link_libraries(FortranHighs PUBLIC highs)
  endif()

  install(TARGETS FortranHighs
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs
    MODULES DESTINATION modules)
  if(NOT MSVC)
    install(FILES ${HIGHS_BINARY_DIR}/modules/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  else()
    install(FILES ${HIGHS_BINARY_DIR}/modules/${CMAKE_BUILD_TYPE}/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  endif()
  # use link rpath
  # set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH
  #   "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif(FORTRAN_FOUND)

if(CSHARP_FOUND)
  message(STATUS "CSharp supported")
  set(csharpsources
    interfaces/highs_csharp_api.cs)
  add_library(HighsCsharp SHARED interfaces/highs_csharp_api.cs)
  add_library(${PROJECT_NAMESPACE}::HighsCsharp ALIAS HighsCsharp)
  target_compile_options(HighsCsharp PUBLIC "/unsafe")
  add_executable(csharpexample ../examples/call_highs_from_csharp.cs)
  target_compile_options(csharpexample PUBLIC "/unsafe")
  target_link_libraries(csharpexample PUBLIC HighsCsharp)
else()
  message(STATUS "No CSharp support")
endif()

find_package(Threads)
if(Threads_FOUND)
  include(CheckAtomic)
  if(HAVE_CXX_ATOMICS64_WITH_LIB)
    if(FAST_BUILD)
      target_link_libraries(highs PRIVATE atomic)
    else()
      target_link_libraries(libhighs atomic)
    endif()
  endif()
endif()

if(FAST_BUILD)
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
      set(THREADS_PREFER_PTHREAD_FLAG ON)
  endif()

  target_link_libraries(highs PRIVATE Threads::Threads)
  target_link_libraries(highs PRIVATE ${CMAKE_DL_LIBS})
else()
  target_link_libraries(libhighs Threads::Threads)
endif()
