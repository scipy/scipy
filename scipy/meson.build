# Platform detection
is_windows = host_machine.system() == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_global_link_arguments(gcc_link_args, language: ['c', 'cpp'])
    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_global_arguments('-mlong-double-64', language: 'c')
    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_global_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('_build_utils/gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_global_arguments('-DMS_WIN64', language: ['c', 'cpp', 'fortran'])
    endif
    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  endif
  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_global_link_arguments(gcc_link_args, language: ['fortran'])
    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_global_arguments('-fno-optimize-sibling-calls',
      language: ['fortran'])
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
  blas_name = ['openblas', 'OpenBLAS']
endif
if lapack_name == 'openblas'
  lapack_name = ['openblas', 'OpenBLAS']
endif
blas = dependency(blas_name)
lapack = dependency(lapack_name)

if blas.name() == 'mkl' or lapack.name() == 'mkl' or get_option('use-g77-abi')
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
  install_dir: py3.get_install_dir() + '/scipy'
)

generate_version = custom_target(
  'generate-version',
  install: true,
  build_always_stale: true,
  build_by_default: true,
  output: 'version.py',
  input: '../tools/version_utils.py',
  command: [py3, '@INPUT@', '--source-root', '@SOURCE_ROOT@'],
  install_dir: py3.get_install_dir() + '/scipy'
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

c_undefined_ok = ['-Wno-maybe-uninitialized']

# Suppress warning for deprecated Numpy API.
# (Suppress warning messages emitted by #warning directives).
# Replace with numpy_nodepr_api after Cython 3.0 is out
cython_c_args += ['-Wno-cpp']
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
