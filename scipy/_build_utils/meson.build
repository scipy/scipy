# These wrappers are only used for _fblas/_flapack/_fblas_64/_flapack_64,
# because the .pyf wrappers for those didn't get source-level symbol rename
# handling (that'd be a lot of work).
#
# These wrappers must be generated gated on the condition that the symbols have
# a prefix and/or a suffix.

blas_wrapper_gen = find_program('_generate_blas_wrapper.py')
needs_lp64_fblas = not use_ilp64 or not cython_blas_ilp64

if (uses_accelerate or blas_name == 'scipy-openblas') and needs_lp64_fblas
  blas_lapack_wrappers_lp64 = custom_target('blas_lapack_wrappers_lp64',
      output: ['blas_lapack_wrappers_lp64.c'],
      command: [blas_wrapper_gen, '-o', '@OUTDIR@', accelerate_flag],
      depend_files: [
        '../linalg/cython_blas_signatures.txt',
        '../linalg/cython_lapack_signatures.txt',
      ]
  )

  blas_lapack_wrapper_lib_lp64 = static_library('blas_lapack_wrapper_lib_lp64',
    blas_lapack_wrappers_lp64,
    c_args: numpy_nodepr_api,
    dependencies: [lapack_lp64_dep, np_dep, py3_dep],
    include_directories: ['../linalg/', './src'],
    gnu_symbol_visibility: 'hidden',
  )
else
  blas_lapack_wrapper_lib_lp64 = []
endif

if use_ilp64 or uses_accelerate or blas_name == 'scipy-openblas'
  blas_lapack_wrappers_ilp64 = custom_target('blas_lapack_wrappers_ilp64',
      output: ['blas_lapack_wrappers_ilp64.c'],
      command: [blas_wrapper_gen, '-o', '@OUTDIR@', accelerate_flag, '--ilp64'],
      depend_files: [
        '../linalg/cython_blas_signatures.txt',
        '../linalg/cython_lapack_signatures.txt',
      ]
  )

  blas_lapack_wrapper_lib_ilp64 = static_library('blas_lapack_wrapper_lib_ilp64',
    blas_lapack_wrappers_ilp64,
    c_args: numpy_nodepr_api,
    dependencies: [lapack_dep, np_dep, py3_dep],
    include_directories: ['../linalg/', './src'],
    gnu_symbol_visibility: 'hidden',
  )
else
  blas_lapack_wrapper_lib_ilp64 = []
endif
