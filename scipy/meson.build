# Platform detection
is_mingw = is_windows and cc.get_id() == 'gcc'
if is_mingw and ff.get_id() != 'gcc'
  error('If you are using GCC on Windows, you must also use GFortran! Detected ' + ff.get_id())
endif

cython_c_args = ['-DCYTHON_CCOMPLEX=0'] # see gh-18975 for why we need this
if is_mingw
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  add_project_link_arguments(gcc_link_args, language: ['c', 'cpp', 'fortran'])
  # Force gcc to float64 long doubles for compatibility with MSVC
  # builds, for C only.
  add_project_arguments('-mlong-double-64', language: 'c')
  # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
  add_project_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
  # Silence warnings emitted by PyOS_snprintf for (%zd), see
  # https://github.com/rgommers/scipy/issues/118.
  # Use as c_args for extensions containing Cython code
  cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  # Flag needed to work around BLAS and LAPACK Gfortran dependence on
  # undocumented C feature when passing single character string arguments. See:
  #   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
  #   https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
  add_project_arguments('-fno-optimize-sibling-calls', language: ['fortran'])
endif

thread_dep = dependency('threads', required: false)

# NumPy include directory - needed in all submodules
# The chdir is needed because within numpy there's an `import signal`
# statement, and we don't want that to pick up scipy's signal module rather
# than the stdlib module. The try-except is needed because when things are
# split across drives on Windows, there is no relative path and an exception
# gets raised. There may be other such cases, so add a catch-all and switch to
# an absolute path. Relative paths are needed when for example a virtualenv is
# placed inside the source tree; Meson rejects absolute paths to places inside
# the source tree.
# For cross-compilation it is often not possible to run the Python interpreter
# in order to retrieve numpy's include directory. It can be specified in the
# cross file instead:
#   [properties]
#   numpy-include-dir = /abspath/to/host-pythons/site-packages/numpy/core/include
#
# This uses the path as is, and avoids running the interpreter.
incdir_numpy = meson.get_external_property('numpy-include-dir', 'not-given')
if incdir_numpy == 'not-given'
  incdir_numpy = run_command(py3,
    [
      '-c',
      '''import os
os.chdir(os.path.join("..", "tools"))
import numpy as np
try:
  incdir = os.path.relpath(np.get_include())
except Exception:
  incdir = np.get_include()
print(incdir)
  '''
    ],
    check: true
  ).stdout().strip()

  # We do need an absolute path to feed to `cc.find_library` below
  _incdir_numpy_abs = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check: true
  ).stdout().strip()
else
  _incdir_numpy_abs = incdir_numpy
endif
inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

npymath_path = _incdir_numpy_abs / '..' / 'lib'
npyrandom_path = _incdir_numpy_abs / '..' / '..' / 'random' / 'lib'
npymath_lib = cc.find_library('npymath', dirs: npymath_path)
npyrandom_lib = cc.find_library('npyrandom', dirs: npyrandom_path)

pybind11_dep = dependency('pybind11', version: '>=2.10.4')

# Pythran include directory and build flags
if use_pythran
  # This external-property may not be needed if we can use the native include
  # dir, see https://github.com/serge-sans-paille/pythran/issues/1394
  incdir_pythran = meson.get_external_property('pythran-include-dir', 'not-given')
  if incdir_pythran == 'not-given'
    incdir_pythran = run_command(py3,
      [
        '-c',
      '''import os
os.chdir(os.path.join("..", "tools"))
import pythran
try:
  incdir = os.path.relpath(pythran.get_include())
except Exception:
  incdir = pythran.get_include()
print(incdir)
'''
      ],
      check: true
    ).stdout().strip()
  endif
  pythran_dep = declare_dependency(
    include_directories: incdir_pythran,
    dependencies: xsimd_dep,
  )
else
  pythran_dep = []
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

# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: numpy_nodepr_api,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

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
  blas = dependency(['openblas', 'OpenBLAS'])
else
  blas = dependency(blas_name)
endif
if blas_name == 'blas'
  # Netlib BLAS has a separate `libcblas.so` which we use directly in the g77
  # ABI wrappers, so detect it and error out if we cannot find it.
  # In the future, this should be done automatically for:
  #   `dependency('blas', modules: cblas)`
  # see https://github.com/mesonbuild/meson/pull/10921.
  cblas = dependency('cblas')
else
  cblas = []
endif

if lapack_name == 'openblas'
  lapack = dependency(['openblas', 'OpenBLAS'])
else
  lapack = dependency(lapack_name)
endif

dependency_map = {
  'BLAS': blas,
  'LAPACK': lapack,
  'PYBIND11': pybind11_dep,
}

# FIXME: conda-forge sets MKL_INTERFACE_LAYER=LP64,GNU, see gh-11812.
#        This needs work on gh-16200 to make MKL robust. We should be
#        requesting `mkl-dynamic-lp64-seq` here. And then there's work needed
#        in general to enable the ILP64 interface (also for OpenBLAS).
uses_mkl = blas.name().to_lower().startswith('mkl') or lapack.name().to_lower().startswith('mkl')
uses_accelerate = blas.name().to_lower().startswith('accelerate') or lapack.name().to_lower().startswith('accelerate')
use_g77_abi = uses_mkl or uses_accelerate or get_option('use-g77-abi')
if use_g77_abi
  g77_abi_wrappers = declare_dependency(
    sources:
      [
        '_build_utils/src/wrap_g77_abi_f.f',
        '_build_utils/src/wrap_g77_abi_c.c'
      ],
    include_directories: inc_np,
    dependencies: [py3_dep, cblas],
  )
else
  g77_abi_wrappers = declare_dependency(sources: ['_build_utils/src/wrap_dummy_g77_abi.f'])
endif

scipy_dir = py3.get_install_dir() / 'scipy'

generate_version = custom_target(
  'generate-version',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: 'version.py',
  input: '../tools/version_utils.py',
  command: [py3, '@INPUT@', '--source-root', '@SOURCE_ROOT@'],
  install_dir: scipy_dir
)

python_sources = [
  '__init__.py',
  '_distributor_init.py',
  'conftest.py',
  'linalg.pxd',
  'optimize.pxd',
  'special.pxd',
]

py3.install_sources(
  python_sources,
  subdir: 'scipy'
)

py3.install_sources(
  ['_build_utils/tests/test_scipy_version.py'],
  subdir: 'scipy/_lib/tests'
)

# Copy the main __init__.py and pxd files to the build dir.
# Needed to trick Cython, it won't do a relative import outside a package
fs = import('fs')
#_cython_tree = declare_dependency(sources: [
_cython_tree = [
  fs.copyfile('__init__.py'),
  fs.copyfile('linalg.pxd'),
  fs.copyfile('optimize.pxd'),
  fs.copyfile('special.pxd'),
]

cython_args = ['-3', '--fast-fail', '--output-file', '@OUTPUT@', '--include-dir', '@BUILD_ROOT@', '@INPUT@']
cython_cplus_args = ['--cplus'] + cython_args

cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : _cython_tree)

cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
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
Wno_incompatible_pointer_types = cc.get_supported_arguments('-Wno-incompatible-pointer-types')

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
_fflag_Wno_surprising = ff.get_supported_arguments('-Wno-surprising')
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

# Intel Fortran (ifort) does not run the preprocessor by default, if Fortran
# code uses preprocessor statements, add this compile flag to it.
_fflag_fpp = []
if ff.get_id() == 'intel-cl'
  if is_windows
    _fflag_fpp = ff.get_supported_arguments('/fpp')
  else
    _fflag_fpp = ff.get_supported_arguments('-fpp')
  endif
endif

# Deal with M_PI & friends; add `use_math_defines` to c_args or cpp_args
# Cython doesn't always get this right itself (see, e.g., gh-16800), so
# explicitly add the define as a compiler flag for Cython-generated code.
if is_windows
  use_math_defines = ['-D_USE_MATH_DEFINES']
else
  use_math_defines = []
endif

# Determine whether it is necessary to link libatomic. This could be the case
# e.g. on 32-bit platforms when atomic operations are used on 64-bit types.
# The check is copied from Mesa <https://www.mesa3d.org/>.
# Note that this dependency is not desired, it came in with a HiGHS update.
# We should try to get rid of it. For discussion, see gh-17777.
null_dep = dependency('', required : false)
atomic_dep = null_dep
code_non_lockfree = '''
  #include <stdint.h>
  int main() {
   struct {
     uint64_t *v;
   } x;
   return (int)__atomic_load_n(x.v, __ATOMIC_ACQUIRE) &
          (int)__atomic_add_fetch(x.v, (uint64_t)1, __ATOMIC_ACQ_REL);
  }
'''
if cc.get_id() != 'msvc'
  if not cc.links(
      code_non_lockfree,
      name : 'Check atomic builtins without -latomic'
    )
    atomic_dep = cc.find_library('atomic', required: false)
    if atomic_dep.found()
      # We're not sure that with `-latomic` things will work for all compilers,
      # so verify and only keep libatomic as a dependency if this works. It is
      # possible the build will fail later otherwise - unclear under what
      # circumstances (compilers, runtimes, etc.) exactly.
      if not cc.links(
          code_non_lockfree,
          dependencies: atomic_dep,
          name : 'Check atomic builtins with -latomic'
        )
        atomic_dep = null_dep
      endif
    endif
  endif
endif

# Suppress warning for deprecated Numpy API.
# (Suppress warning messages emitted by #warning directives).
# Replace with numpy_nodepr_api after Cython 3.0 is out
cython_c_args += [_cpp_Wno_cpp, use_math_defines]
cython_cpp_args = cython_c_args

compilers = {
  'C': cc,
  'CPP': cpp,
  'CYTHON': meson.get_compiler('cython'),
  'FORTRAN': meson.get_compiler('fortran')
}

machines = {
  'HOST': host_machine,
  'BUILD': build_machine,
}

conf_data = configuration_data()

# Set compiler information
foreach name, compiler : compilers
  conf_data.set(name + '_COMP', compiler.get_id())
  conf_data.set(name + '_COMP_LINKER_ID', compiler.get_linker_id())
  conf_data.set(name + '_COMP_VERSION', compiler.version())
  conf_data.set(name + '_COMP_CMD_ARRAY', ', '.join(compiler.cmd_array()))
  conf_data.set(name + '_COMP_ARGS', ', '.join(
      get_option(name.to_lower() + '_args')
    )
  )
  conf_data.set(name + '_COMP_LINK_ARGS', ', '.join(
      get_option(name.to_lower() + '_link_args')
    )
  )
endforeach
# Add `pythran` information if present
if use_pythran
  conf_data.set('PYTHRAN_VERSION', pythran.version())
  conf_data.set('PYTHRAN_INCDIR', incdir_pythran)
endif

# Machines CPU and system information
foreach name, machine : machines
  conf_data.set(name + '_CPU', machine.cpu())
  conf_data.set(name + '_CPU_FAMILY', machine.cpu_family())
  conf_data.set(name + '_CPU_ENDIAN', machine.endian())
  conf_data.set(name + '_CPU_SYSTEM', machine.system())
endforeach

conf_data.set('CROSS_COMPILED', meson.is_cross_build())

# Python information
conf_data.set('PYTHON_PATH', py3.full_path())
conf_data.set('PYTHON_VERSION', py3.language_version())

# Dependencies information
foreach name, dep : dependency_map
  conf_data.set(name + '_NAME', dep.name())
  conf_data.set(name + '_FOUND', dep.found())
  if dep.found()
    conf_data.set(name + '_VERSION', dep.version())
    conf_data.set(name + '_TYPE_NAME', dep.type_name())
    conf_data.set(name + '_INCLUDEDIR', dep.get_variable('includedir', default_value: 'unknown'))
    conf_data.set(name + '_LIBDIR', dep.get_variable('libdir', default_value: 'unknown'))
    conf_data.set(name + '_OPENBLAS_CONFIG', dep.get_variable('openblas_config', default_value: 'unknown'))
    conf_data.set(name + '_PCFILEDIR', dep.get_variable('pcfiledir', default_value: 'unknown'))
  endif
endforeach

configure_file(
  input: '__config__.py.in',
  output: '__config__.py',
  configuration : conf_data,
  install_dir: scipy_dir,
)

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
subdir('io')
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
subdir('datasets')
subdir('misc')
