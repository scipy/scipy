if generate_blas_wrappers
  blas_wrapper_gen = find_program('_generate_blas_wrapper.py')
  blas_lapack_wrappers = custom_target('blas_lapack_wrappers',
      output: ['blas_lapack_wrappers.c'],
      command: [blas_wrapper_gen, '-o', '@OUTDIR@', accelerate_flag],
      depend_files: [
        '../linalg/cython_blas_signatures.txt',
        '../linalg/cython_lapack_signatures.txt',
      ]
  )

  blas_lapack_wrapper_lib = static_library('blas_lapack_wrapper_lib',
    blas_lapack_wrappers,
    c_args: numpy_nodepr_api,
    dependencies: [blas_dep, lapack_dep, np_dep, py3_dep],
    include_directories: ['../linalg/', './src'],
    gnu_symbol_visibility: 'hidden',
  )
else
  blas_lapack_wrapper_lib = []
endif
