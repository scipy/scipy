# Platform detection
is_mingw = is_windows and cc.get_id() == 'gcc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_project_link_arguments(gcc_link_args, language: ['c', 'cpp'])
    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_project_arguments('-mlong-double-64', language: 'c')
    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_project_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('_build_utils/gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_project_arguments('-DMS_WIN64', language: ['c', 'cpp', 'fortran'])
    endif
    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  endif
  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_project_link_arguments(gcc_link_args, language: ['fortran'])
    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_project_arguments('-fno-optimize-sibling-calls', language: ['fortran'])
  endif
endif

# NumPy include directory - needed in all submodules
incdir_numpy = run_command(py3,
  [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

cc = meson.get_compiler('c')
npymath_path = incdir_numpy / '..' / 'lib'
npymath_lib = cc.find_library('npymath', dirs: npymath_path)
npyrandom_path = incdir_numpy / '..' / '..' / 'random' / 'lib'
# Note: `required: false` can be removed once numpy 1.19 is the minimum version
npyrandom_lib = cc.find_library('npyrandom', dirs: npyrandom_path, required: false)

# pybind11 include directory - needed in several submodules
incdir_pybind11 = run_command(py3,
  [
    '-c',
    'import pybind11; print(pybind11.get_include())'
  ],
  check: true
).stdout().strip()

inc_pybind11 = include_directories(incdir_pybind11)

# Pythran include directory and build flags
use_pythran = run_command(py3,
  [
    '-c',
    'import os; print(os.environ.get("SCIPY_USE_PYTHRAN", 1))'
  ],
  check: true
).stdout().strip() == '1'

if use_pythran
  incdir_pythran = run_command(py3,
    [
      '-c',
      'import os; os.chdir(".."); import pythran; print(os.path.dirname(pythran.__file__));'
    ],
    check: true
  ).stdout().strip()
  inc_pythran = include_directories(incdir_pythran)
else
  inc_pythran = []
endif

# Note: warning flags are added to this further down
cpp_args_pythran = [
  '-DENABLE_PYTHON_MODULE',
  '-D__PYTHRAN__=3',
  '-DPYTHRAN_BLAS_NONE'
]

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'


# TODO: 64-bit BLAS and LAPACK
#
# Note that this works as long as BLAS and LAPACK are detected properly via
# pkg-config. By default we look for OpenBLAS, other libraries can be configured via
# `meson configure -Dblas=blas -Dlapack=lapack` (example to build with Netlib
# BLAS and LAPACK).
# For MKL and for auto-detecting one of multiple libs, we'll need a custom
# dependency in Meson (like is done for scalapack) - see
# https://github.com/mesonbuild/meson/issues/2835
blas_name = get_option('blas')
lapack_name = get_option('lapack')
# pkg-config uses a lower-case name while CMake uses a capitalized name, so try
# that too to make the fallback detection with CMake work
if blas_name == 'openblas'
  blas_name = ['openblas', 'OpenBLAS']
endif
if lapack_name == 'openblas'
  lapack_name = ['openblas', 'OpenBLAS']
endif
blas = dependency(blas_name)
lapack = dependency(lapack_name)

# FIXME: conda-forge sets MKL_INTERFACE_LAYER=LP64,GNU, see gh-11812.
#        This needs work on gh-16200 to make MKL robust. We should be
#        requesting `mkl-dynamic-lp64-seq` here. And then there's work needed
#        in general to enable the ILP64 interface (also for OpenBLAS).
uses_mkl = blas.name().to_lower().startswith('mkl') or lapack.name().to_lower().startswith('mkl')
if uses_mkl or get_option('use-g77-abi')
  g77_abi_wrappers = files([
    '_build_utils/src/wrap_g77_abi_f.f',
    '_build_utils/src/wrap_g77_abi_c.c'
    ])
else
  g77_abi_wrappers = files('_build_utils/src/wrap_dummy_g77_abi.f')
endif

generate_config = custom_target(
  'generate-config',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: '__config__.py',
  input: '../tools/config_utils.py',
  command: [py3, '@INPUT@', '@OUTPUT@'],
  install_dir: py3.get_install_dir(pure: false) / 'scipy'
)

generate_version = custom_target(
  'generate-version',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: 'version.py',
  input: '../tools/version_utils.py',
  command: [py3, '@INPUT@', '--source-root', '@SOURCE_ROOT@'],
  install_dir: py3.get_install_dir(pure: false) / 'scipy'
)

python_sources = [
  '__init__.py',
  '_distributor_init.py',
  'conftest.py',
  'linalg.pxd',
  'optimize.pxd',
  'special.pxd'
]

py3.install_sources(
  python_sources,
  pure: false,
  subdir: 'scipy'
)

py3.install_sources(
  ['_build_utils/tests/test_scipy_version.py'],
  pure: false,
  subdir: 'scipy/_lib/tests'
)

# Needed to trick Cython, it won't do a relative import outside a package
_cython_tree = custom_target('_cython_tree',
  output: [
    '__init__.py',
    'linalg.pxd',
    'optimize.pxd',
    'special.pxd'
  ],
  input: [
    '__init__.py',
    'linalg.pxd',
    'optimize.pxd',
    'special.pxd'
  ],
  command: [copier, '@INPUT@', '@OUTDIR@'],
)
cython_tree = declare_dependency(sources: _cython_tree)

cython_cli = find_program('_build_utils/cythoner.py')

cython_gen = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@'],
  output : '@BASENAME@.c',
  depends : _cython_tree)

cython_gen_cpp = generator(cython_cli,
  arguments : ['@INPUT@', '@OUTPUT@', '--cplus'],
  output : '@BASENAME@.cpp',
  depends : [_cython_tree])

# Check if compiler flags are supported. This is necessary to ensure that SciPy
# can be built with any supported compiler. We need so many warning flags
# because we want to be able to build with `-Werror` in CI; that ensures that
# for new code we add, there are no unexpected new issues introduced.
#
# Cleaning up code so we no longer need some of these warning flags is useful,
# but not a priority.
#
# The standard convention used here is:
#   - for C, drop the leading dash and turn remaining dashes into underscores
#   - for C++, prepend `_cpp` and turn remaining dashes into underscores
#   - for Fortran, prepend `_fflags` and turn remaining dashes into underscores

# C warning flags
Wno_maybe_uninitialized = cc.get_supported_arguments('-Wno-maybe-uninitialized')
Wno_discarded_qualifiers = cc.get_supported_arguments('-Wno-discarded-qualifiers')
Wno_empty_body = cc.get_supported_arguments('-Wno-empty-body')
Wno_implicit_function_declaration = cc.get_supported_arguments('-Wno-implicit-function-declaration')
Wno_parentheses = cc.get_supported_arguments('-Wno-parentheses')
Wno_switch = cc.get_supported_arguments('-Wno-switch')
Wno_unused_label = cc.get_supported_arguments('-Wno-unused-label')
Wno_unused_variable = cc.get_supported_arguments('-Wno-unused-variable')

# C++ warning flags
_cpp_Wno_cpp = cpp.get_supported_arguments('-Wno-cpp')
_cpp_Wno_deprecated_declarations = cpp.get_supported_arguments('-Wno-deprecated-declarations')
_cpp_Wno_class_memaccess = cpp.get_supported_arguments('-Wno-class-memaccess')
_cpp_Wno_format_truncation = cpp.get_supported_arguments('-Wno-format-truncation')
_cpp_Wno_non_virtual_dtor = cpp.get_supported_arguments('-Wno-non-virtual-dtor')
_cpp_Wno_sign_compare = cpp.get_supported_arguments('-Wno-sign-compare')
_cpp_Wno_switch = cpp.get_supported_arguments('-Wno-switch')
_cpp_Wno_terminate = cpp.get_supported_arguments('-Wno-terminate')
_cpp_Wno_unused_but_set_variable = cpp.get_supported_arguments('-Wno-unused-but-set-variable')
_cpp_Wno_unused_function = cpp.get_supported_arguments('-Wno-unused-function')
_cpp_Wno_unused_local_typedefs = cpp.get_supported_arguments('-Wno-unused-local-typedefs')
_cpp_Wno_unused_variable = cpp.get_supported_arguments('-Wno-unused-variable')
_cpp_Wno_int_in_bool_context = cpp.get_supported_arguments('-Wno-int-in-bool-context')

cpp_args_pythran += [
  _cpp_Wno_cpp,
  _cpp_Wno_deprecated_declarations,
  _cpp_Wno_unused_but_set_variable,
  _cpp_Wno_unused_function,
  _cpp_Wno_unused_variable,
  _cpp_Wno_int_in_bool_context,
]

# Fortran warning flags
_fflag_Wno_argument_mismatch = ff.get_supported_arguments('-Wno-argument-mismatch')
_fflag_Wno_conversion = ff.get_supported_arguments('-Wno-conversion')
_fflag_Wno_intrinsic_shadow = ff.get_supported_arguments('-Wno-intrinsic-shadow')
_fflag_Wno_maybe_uninitialized = ff.get_supported_arguments('-Wno-maybe-uninitialized')
_fflag_Wno_uninitialized = ff.get_supported_arguments('-Wno-uninitialized')
_fflag_Wno_unused_dummy_argument = ff.get_supported_arguments('-Wno-unused-dummy-argument')
_fflag_Wno_unused_label = ff.get_supported_arguments('-Wno-unused-label')
_fflag_Wno_unused_variable = ff.get_supported_arguments('-Wno-unused-variable')
_fflag_Wno_tabs = ff.get_supported_arguments('-Wno-tabs')
# The default list of warnings to ignore from Fortran code. There is a lot of
# old, vendored code that is very bad and we want to compile it silently (at
# least with GCC and Clang)
fortran_ignore_warnings = ff.get_supported_arguments(
 _fflag_Wno_argument_mismatch,
 _fflag_Wno_conversion,
 _fflag_Wno_maybe_uninitialized,
 _fflag_Wno_unused_dummy_argument,
 _fflag_Wno_unused_label,
 _fflag_Wno_unused_variable,
 _fflag_Wno_tabs,
)

# Deal with M_PI & friends; add `use_math_defines` to c_args or cpp_args
# Cython doesn't always get this right itself (see, e.g., gh-16800), so
# explicitly add the define as a compiler flag for Cython-generated code.
if is_windows
  use_math_defines = ['-D_USE_MATH_DEFINES']
else
  use_math_defines = []
endif

# Suppress warning for deprecated Numpy API.
# (Suppress warning messages emitted by #warning directives).
# Replace with numpy_nodepr_api after Cython 3.0 is out
cython_c_args += ['-Wno-cpp', use_math_defines]
cython_cpp_args = cython_c_args

# Ordering of subdirs: special and linalg come first, because other submodules
# have dependencies on cython_special.pxd and cython_linalg.pxd. After those,
# subdirs with the most heavy builds should come first (that parallelizes
# better)
subdir('_lib')
subdir('special')
subdir('linalg')
subdir('sparse')
subdir('stats')
subdir('fft')
subdir('spatial')
subdir('cluster')
subdir('constants')
subdir('fftpack')
subdir('integrate')
subdir('signal')
subdir('interpolate')
subdir('ndimage')
subdir('odr')
subdir('optimize')
subdir('misc')
subdir('io')
