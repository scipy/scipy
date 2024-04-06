cpp_args = [
  '-DBOOST_MATH_STANDALONE=1',
  '-DBOOST_MATH_DOMAIN_ERROR_POLICY=ignore_error',
  '-DBOOST_MATH_EVALUATION_ERROR_POLICY=user_error',
  '-DBOOST_MATH_OVERFLOW_ERROR_POLICY=user_error',
  '-DBOOST_MATH_PROMOTE_DOUBLE_POLICY=false'
]

pyx_files = [
  ['binom_ufunc', binom_ufunc_pyx],
  ['beta_ufunc', beta_ufunc_pyx],
  ['hypergeom_ufunc', hypergeom_ufunc_pyx],
  ['nbinom_ufunc', nbinom_ufunc_pyx],
  ['ncf_ufunc', ncf_ufunc_pyx],
  ['ncx2_ufunc', ncx2_ufunc_pyx],
  ['nct_ufunc', nct_ufunc_pyx],
  ['skewnorm_ufunc', skewnorm_ufunc_pyx],
  ['invgauss_ufunc', invgauss_ufunc_pyx],
]

foreach pyx_file: pyx_files
  py3.extension_module(pyx_file[0],
    pyx_file[1],
    include_directories: ['include', '../../_lib/boost_math/include'],
    dependencies: np_dep,
    link_args: version_link_args,
    cpp_args: [cpp_args, cython_cpp_args],
    install: true,
    subdir: 'scipy/stats/_boost'
  )
endforeach

py3.install_sources([
    '__init__.py'
  ],
  subdir: 'scipy/stats/_boost'
)
