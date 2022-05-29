project(
  'SciPy',
  'c', 'cpp',
  # Note that the git commit hash cannot be added dynamically here (it is added
  # in the dynamically generated and installed `scipy/version.py` though - see
  # tools/version_utils.py
  version: '1.9.0.dev0',
  license: 'BSD-3',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'cpp_std=c++14',  # TODO: use c++11 if 14 is not available
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',
    'fortran_std=legacy',
    'blas=openblas',
    'lapack=openblas'
  ],
)

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
# This argument is called -Wno-unused-but-set-variable by GCC, however Clang
# doesn't recognize that.
if cc.has_argument('-Wno-unused-but-set-variable')
  add_global_arguments('-Wno-unused-but-set-variable', language : 'c')
endif

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

cython = find_program('cython')
pythran = find_program('pythran')
generate_f2pymod = files('tools/generate_f2pymod.py')
tempita = files('scipy/_build_utils/tempita.py')
copier = find_program(['cp', 'scipy/_build_utils/copyfiles.py'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()

subdir('scipy')
