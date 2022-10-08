project(
  'SciPy',
  'c', 'cpp',
  # Note that the git commit hash cannot be added dynamically here (it is added
  # in the dynamically generated and installed `scipy/version.py` though - see
  # tools/version_utils.py
  version: '1.9.2',
  license: 'BSD-3',
  meson_version: '>= 0.62.2',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'cpp_std=c++14',
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

# TODO: the below -Wno flags are all needed to silence warnings in
# f2py-generated code. This should be fixed in f2py itself.
_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
  '-Wno-incompatible-pointer-types',
)
add_project_arguments(_global_c_args, language : 'c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

is_windows = host_machine.system() == 'windows'

# Intel Fortran on Windows does things differently, so deal with that
if is_windows and ff.get_id() == 'intel-cl'
  _ifort_flags = ff.get_supported_arguments('/MD', '/names:lowercase', '/assume:underscore')
  add_project_arguments(_ifort_flags, language: 'fortran')
endif

cython = find_program('cython')
pythran = find_program('pythran')
generate_f2pymod = files('tools/generate_f2pymod.py')
tempita = files('scipy/_build_utils/tempita.py')
copier = find_program(['cp', 'scipy/_build_utils/copyfiles.py'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
# NOTE: with Meson >=0.64.0 we can add `pure: false` here and remove that line
# everywhere else, see https://github.com/mesonbuild/meson/pull/10783.
py3 = py_mod.find_installation()
# SciPy 1.9.0-specific error message, see the same message in
# `scipy/__init__.py` and gh-14986
if py3.language_version().version_compare('>=3.12')
  error('Your Python version is too new. SciPy 1.9 supports ' +
        'Python 3.8-3.11; if you are trying to build from source for the ' +
        'most recent SciPy version you may hit this error as well. Please ' +
        'build from the `main` branch on GitHub instead.')
endif
py3_dep = py3.dependency()

subdir('scipy')
