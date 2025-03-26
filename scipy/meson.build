fs = import('fs')

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


# Uses the `numpy-config` executable (or a user's numpy.pc pkg-config file).
# Will work for numpy>=2.0, hence not required (it'll be a while until 2.0 is
# our minimum supported version). Using this now to be able to detect the
# version easily for >=2.0.
_numpy_dep = dependency('numpy', required: false)
f2py_freethreading_arg = []
if _numpy_dep.found()
  if _numpy_dep.version().version_compare('>=2.1.0')
    f2py_freethreading_arg = ['--free-threading']
  endif
endif

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
# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation.
numpy_nodepr_api = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']
np_dep = declare_dependency(include_directories: inc_np, compile_args: numpy_nodepr_api)


incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

npymath_path = _incdir_numpy_abs / '..' / 'lib'
npyrandom_path = _incdir_numpy_abs / '..' / '..' / 'random' / 'lib'
npymath_lib = cc.find_library('npymath', dirs: npymath_path)
npyrandom_lib = cc.find_library('npyrandom', dirs: npyrandom_path)

pybind11_dep = dependency('pybind11', version: '>=2.13.2')

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

# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
  gnu_symbol_visibility: 'hidden',
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

cdata = configuration_data()

# Test variable attribute to use for thread-local storage;
# Adapted from `numpy/_core/meson.build`.
check_tls_attrs = [
  ['thread_local', 'HAVE_THREAD_LOCAL'],    # C23
  ['_Thread_local', 'HAVE__THREAD_LOCAL'],  # C11/C17
  ['__thread', 'HAVE__THREAD'],
]
if is_windows and not is_mingw
  check_tls_attrs += ['__declspec(thread)', 'HAVE___DECLSPEC_THREAD_']
endif
f2py_tls_define = ''
foreach tls_attrs: check_tls_attrs
  attr = tls_attrs[0]
  code = f'''
    #pragma GCC diagnostic error "-Wattributes"
    #pragma clang diagnostic error "-Wattributes"

    int @attr@ foo;
  '''
  code += '''
    int
    main()
    {
      return 0;
    }
  '''
  if cc.compiles(code, name: tls_attrs[0])
    cdata.set10(tls_attrs[1], true)
    f2py_tls_define = tls_attrs[0]
  endif
endforeach

# Contains only TLS check results for now - name chosen for when more compiler
# checks need adding in the future.
scipy_config_h = configure_file(
  input: 'scipy_config.h.in',
  output: 'scipy_config.h',
  configuration: cdata,
  install: false
)

_f2py_c_args = [f'-DF2PY_THREAD_LOCAL_DECL=@f2py_tls_define@']
fortranobject_dep = declare_dependency(
  dependencies: fortranobject_dep,
  compile_args: _f2py_c_args,
)

f2py = find_program('f2py')
# It should be quite rare for the `f2py` executable to not be the one from
# `numpy` installed in the Python env we are building for (unless we are
# cross-compiling). If it is from a different env, that is still fine as long
# as it's not too old. We are only using f2py as a code generator, and the
# output is not dependent on platform or Python version (see gh-20612 for more
# details).
# This should be robust enough. If not, we can make this more complex, using
# a fallback to `python -m f2py` rather than erroring out.
f2py_version = run_command([f2py, '-v'], check: true).stdout().strip()
if f2py_version.version_compare('<'+min_numpy_version)
  error(f'Found f2py executable is too old: @f2py_version@')
endif

# Note: this generator cannot handle:
# 1. `.pyf.src` files, because `@BASENAME@` will still include .pyf
# 2. targets with #include's (due to no `depend_files` - see feature request
#    at meson#8295)
f2py_gen = generator(generate_f2pymod,
  arguments : ['@INPUT@', '-o', '@BUILD_DIR@', '@EXTRA_ARGS@'] + f2py_freethreading_arg,
  output : ['_@BASENAME@module.c', '_@BASENAME@-f2pywrappers.f'],
)


# Start of BLAS/LAPACK detection

blas_name = get_option('blas')
lapack_name = get_option('lapack')
blas_symbol_suffix = get_option('blas-symbol-suffix')
use_ilp64 = get_option('use-ilp64')
blas_order = get_option('blas-order')
lapack_order = get_option('lapack-order')

# Detect LP64 (32-bit) BLAS and LAPACK, we always need that
blas_interface = ['interface: lp64']
if blas_order == ['auto']
  blas_order = []
  if host_machine.system() == 'darwin'
    blas_order += 'accelerate'
  endif
  if host_machine.cpu_family() == 'x86_64'
    blas_order += 'mkl'
  endif
  blas_order += ['openblas', 'flexiblas', 'blis', 'blas']
endif
if lapack_order == ['auto']
  lapack_order = []
  if host_machine.system() == 'darwin'
    lapack_order += 'accelerate'
  endif
  if host_machine.cpu_family() == 'x86_64'
    lapack_order += 'mkl'
  endif
  lapack_order += ['openblas', 'flexiblas', 'lapack']
endif

# MKL-specific options
_threading_opt = get_option('mkl-threading')
if _threading_opt == 'auto'
  # Switch default to iomp once conda-forge missing openmp.pc issue is fixed
  mkl_opts = ['threading: seq']
else
  mkl_opts = ['threading: ' + _threading_opt]
endif
blas_opts = {'mkl': mkl_opts}
mkl_version_req = '>=2023.0'  # see gh-24824
mkl_may_use_sdl = not use_ilp64 and _threading_opt in ['auto', 'iomp']


# BLAS detection.
#
# First try scipy-openblas, and if found don't look for blas/lapack, we
# know what's inside the scipy-openblas wheels already.
if blas_name == 'openblas' or blas_name == 'auto'
  blas = dependency('scipy-openblas', method: 'pkg-config', required: false)
  if blas.found()
    blas_name = 'scipy-openblas'
  endif
endif
if blas_name == 'auto'
  foreach _name : blas_order
    if _name == 'mkl'
      blas = dependency('mkl',
        modules: blas_interface + mkl_opts,
        required: false,  # may be required, but we need to emit a custom error message
        version: mkl_version_req,
      )
      # Insert a second try with MKL, because we may be rejecting older versions
      # or missing it because no pkg-config installed. If so, we need to retry
      # with MKL SDL, and drop the version constraint (this always worked).
      if not blas.found() and mkl_may_use_sdl
        blas = dependency('mkl', modules: ['sdl: true'], required: false)
      endif
    else
      blas = dependency(_name, modules: blas_interface, required: false)
    endif
    if blas.found()
      break
    endif
  endforeach
else
  if blas_name.to_lower().startswith('mkl')
    mkl_uses_sdl = false
    blas = dependency('mkl',
      modules: blas_interface + mkl_opts,
      required: false,
      version: mkl_version_req,
    )
    # Same deal as above - try again for MKL
    if not blas.found() and mkl_may_use_sdl
      mkl_uses_sdl = true
      blas = dependency('mkl', modules: ['sdl: true'], required: false)
    endif
  else
    blas = dependency(blas_name, modules: blas_interface, required: false)
  endif
endif

if not blas.found()
  error('No BLAS library detected! SciPy needs one, please install it.')
endif

_args_blas = []  # note: used for C and C++ via `blas_dep` below
if blas_symbol_suffix == 'auto'
  blas_symbol_suffix = blas.get_variable('symbol_suffix', default_value: '')
  message(f'BLAS symbol suffix: @blas_symbol_suffix@')
endif
_blas_incdir = []
if blas_symbol_suffix != ''
  _args_blas += ['-DBLAS_SYMBOL_SUFFIX=' + blas_symbol_suffix]
  _blas_incdir = ['.']  # for the generated symbol rename shims in the build dir
endif
blas_dep = declare_dependency(
  dependencies: [blas],
  compile_args: _args_blas,
  include_directories: _blas_incdir,
)

# LAPACK detection
if 'mkl' in blas.name() or blas.name() == 'accelerate' or blas_name == 'scipy-openblas'
  # For these libraries we know that they contain LAPACK, and it's desirable to
  # use that - no need to run the full detection twice.
  lapack = blas
else
  if lapack_name == 'auto'
    foreach _name : lapack_order
      lapack = dependency(_name, modules: ['lapack'] + blas_interface, required: false)
      if lapack.found()
        break
      endif
    endforeach
  else
    lapack = dependency(lapack_name, modules: ['lapack'] + blas_interface, required: false)
  endif
endif

if not lapack.found()
  error('No LAPACK library detected! SciPy needs one, please install it.')
endif
lapack_dep = declare_dependency(dependencies: [lapack, blas_dep])



macOS13_3_or_later = false
if host_machine.system() == 'darwin'
  r = run_command('xcrun', '-sdk', 'macosx', '--show-sdk-version', check: true)
  sdkVersion = r.stdout().strip()

  macOS13_3_or_later = sdkVersion.version_compare('>=13.3')
endif

_args_blas_lapack = []
accelerate_flag = []
generate_blas_wrappers = false
if blas_name == 'accelerate'
  if not macOS13_3_or_later
    error('macOS Accelerate is only supported on macOS >=13.3')
  endif
  if cc.get_id() != 'clang'
    warning('accelerate may not be properly detected with non-native Apple compiler due to https://github.com/mesonbuild/meson/issues/13608')
  endif
  _args_blas_lapack += ['-DACCELERATE_NEW_LAPACK']
  generate_blas_wrappers = true
  accelerate_flag = '-a'
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


dependency_map = {
  'BLAS': blas,
  'LAPACK': lapack,
  'PYBIND11': pybind11_dep,
}

# NB: from this point on blas_name is e.g. 'mkl-lp64-dynamic-seq'
blas_name = blas.name()
lapack_name = lapack.name()
uses_mkl = blas_name.to_lower().startswith('mkl')
uses_accelerate = blas_name.to_lower().startswith('accelerate')
use_g77_abi = uses_mkl or uses_accelerate or get_option('use-g77-abi')
if use_g77_abi
  g77_abi_wrappers = static_library(
    'g77_abi_wrappers',
    ['_build_utils/src/wrap_g77_abi.c'],
    dependencies: [py3_dep, blas, np_dep],
    gnu_symbol_visibility: 'hidden',
  )
else
  g77_abi_wrappers = static_library(
    'dummy_g77_abi_wrappers',
    ['_build_utils/src/wrap_dummy_g77_abi.c'],
    dependencies: [py3_dep, blas, np_dep],
    gnu_symbol_visibility: 'hidden',
  )
endif

# Reuse the names, so we ensure we don't lose the arguments wrapped in with
# declare_dependency. Also, avoids changing `dependencies: blas` to blas_dep in other files.
blas = blas_dep
lapack = lapack_dep

if use_ilp64
  # Okay, we need ILP64 BLAS and LAPACK *in addition to LP64*. So we need to
  # detect the ILP64 variants of the found LP64 libraries now.
  _args_blas = []
  blas_interface = ['interface: ilp64']
  if blas_name == 'flexiblas'
    blas_name = 'flexiblas64'
  endif
  _args_blas += ['-DHAVE_BLAS_ILP64']
  if 'openblas' in blas_name
    _args_blas += ['-DOPENBLAS_ILP64_NAMING_SCHEME']
  endif

  # Run the detection
  if uses_mkl
    mkl_uses_sdl = false   # FIXME, why
    if mkl_uses_sdl
      mkl_opts = ['sdl: true']
    endif
    blas_ilp64 = dependency('mkl', modules: ['interface: ilp64'] + mkl_opts)
    lapack_ilp64 = blas_ilp64
  else
    if blas_name == 'openblas'
      # We cannot allow plain openblas here, that's already the LP64 library and
      # will lead to problems (there is, as of now, no combined OpenBLAS build
      # with 32 and 64 bit symbols)
      blas_name = ['openblas64', 'openblas_ilp64']
    endif
    blas_ilp64 = dependency(blas_name, modules: blas_interface)
    lapack_ilp64 = dependency(lapack_name, modules: ['lapack'] + blas_interface)
  endif

  # Pick up the symbol suffix, it may be auto-detected by Meson and different from LP64
  if blas_symbol_suffix == 'auto'
    if blas_name == 'scipy-openblas'
      blas_symbol_suffix = '64_'
    else
      blas_symbol_suffix = blas_ilp64.get_variable('symbol_suffix', default_value: '')
    endif
    message(f'BLAS symbol suffix (ILP64): @blas_symbol_suffix@')
  endif
  _blas_incdir = []
  if blas_symbol_suffix != ''
    _args_blas += ['-DBLAS_SYMBOL_SUFFIX=' + blas_symbol_suffix]
    _blas_incdir = ['.']
  endif
  # Declare the ILP64 dependencies
  message('BLAS / LAPACK ILP64 detected: ', blas_ilp64.name(), ', ', lapack_ilp64.name())

  blas_ilp64 = declare_dependency(
    dependencies: [blas_ilp64],
    compile_args: _args_blas,
    include_directories: _blas_incdir,
  )
  lapack_ilp64 = declare_dependency(dependencies: [lapack_ilp64, blas_ilp64])
else
  # we're not using ILP64; user will link to the always-available LP64 blas/lapack
  # (all users must use preprocessor macros BLAS_NAME to handle the two options)
  blas_ilp64 = blas
  lapack_ilp64 = lapack
  message('LAPACK ILP64 not requested / not detected.')
endif

# End of BLAS/LAPACK handling


scipy_dir = py3.get_install_dir() / 'scipy'

# Generate version.py for sdist
meson.add_dist_script(
   ['../tools/gitversion.py', '--meson-dist', '--write', 'scipy/version.py']
)
if not fs.exists('version.py')
  generate_version = custom_target(
    'generate-version',
    install: true,
    build_always_stale: true,
    build_by_default: true,
    output: 'version.py',
    input: '../tools/gitversion.py',
    command: ['../tools/gitversion.py', '--meson-dist', '--write', 'scipy/version.py'],
    install_dir: scipy_dir,
    install_tag: 'python-runtime',
  )
else
  # When building from sdist, version.py exists and should be included
  py3.install_sources(['version.py'], subdir : 'scipy')
endif

python_sources = [
  '__init__.py',
  '_distributor_init.py',
  'conftest.py',
 # 'linalg.pxd',
 # 'optimize.pxd',
 # 'special.pxd',
]

if fs.exists('_distributor_init_local.py')
  python_sources += ['_distributor_init_local.py']
endif

py3.install_sources(
  python_sources,
  subdir: 'scipy'
)

# Copy the main __init__.py to the build dir.
# Some submodules (linalg, special, optimize) add pxd files to this.
# Needed to trick Cython, it won't do a relative import outside a package
#_cython_tree = declare_dependency(sources: [
_cython_tree = [fs.copyfile('__init__.py')]

cython_args = ['-3', '--fast-fail', '--output-file', '@OUTPUT@', '--include-dir', '@BUILD_ROOT@', '@INPUT@']
if cy.version().version_compare('>=3.1.0')
  cython_args += ['-Xfreethreading_compatible=True']
endif
cython_cplus_args = ['--cplus'] + cython_args

cython_gen = generator(cython,
  arguments : cython_args,
  output : '@BASENAME@.c',
  depends : _cython_tree)

cython_gen_cpp = generator(cython,
  arguments : cython_cplus_args,
  output : '@BASENAME@.cpp',
  depends : [_cython_tree])

if use_pythran
  # TODO: add argument to mark extension modules as safe to run without the GIL,
  #       once pythran supports that.
  pythran_gen = generator(pythran,
    arguments : ['-E', '@INPUT@', '-o', '@OUTPUT@'],
    output : '@BASENAME@.cpp',
  )
endif

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
Wno_unused_result = cc.get_supported_arguments('-Wno-unused-result')
Wno_unused_variable = cc.get_supported_arguments('-Wno-unused-variable')
Wno_unused_but_set_variable = cc.get_supported_arguments('-Wno-unused-but-set-variable')

# C++ warning flags
_cpp_Wno_bitwise_instead_of_logical = cpp.get_supported_arguments('-Wno-bitwise-instead-of-logical')
_cpp_Wno_cpp = cpp.get_supported_arguments('-Wno-cpp')
_cpp_Wno_class_memaccess = cpp.get_supported_arguments('-Wno-class-memaccess')
_cpp_Wno_deprecated_declarations = cpp.get_supported_arguments('-Wno-deprecated-declarations')
_cpp_Wno_deprecated_builtins = cpp.get_supported_arguments('-Wno-deprecated-builtins')
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

# Gfortran does run the preprocessor for .F files, and PROPACK is the only
# component which needs the preprocessor (unless we need symbol renaming for
# blas_symbol_suffix).
_fflag_preprocess = []
_gfortran_preprocess = ['-cpp', '-ffree-line-length-none', '-ffixed-line-length-none']
if ff.has_multi_arguments(_gfortran_preprocess)
  _fflag_preprocess = _gfortran_preprocess
else
  _fflag_preprocess = ff.first_supported_argument(['-fpp', '/fpp', 'cpp'])
endif

_fflag_lp64 = []
_fflag_ilp64 = []
f2py_ilp64_opts = []
if use_ilp64
  # Gfortran and Clang use `-fdefault-integer-8` to switch to 64-bit integers by
  # default, all other known compilers use `-i8`
  _fflag_ilp64 = ff.first_supported_argument(['-fdefault-integer-8', '-i8'])

  # Write out a mapping file for f2py for defaulting to ILP64
  conf_data = configuration_data()
  if cc.sizeof('long') == 8
    conf_data.set('int64_name', 'long')
  elif cc.sizeof('long long') == 8
    conf_data.set('int64_name', 'long long')
  else
    error('Neither `long` nor `long long` is 64-bit, giving up.')
  endif
  int64_f2cmap = configure_file(
    input: '_build_utils/int64.f2cmap.in',
    output: 'int64.f2cmap',
    configuration: conf_data,
    install: false,
  )
  f2py_ilp64_opts = ['--f2cmap', int64_f2cmap]
endif

if blas_symbol_suffix != ''
  # We need to patch source files that use BLAS/LAPACK symbols.
  # In addition, we now need to enable the Fortran preprocessor on all targets
  # that depend on BLAS/LAPACK
  # Note: this came from `scipy/build_utils/_fortran.py` in the distutils build.
  _fflag_lp64 += _fflag_preprocess  # TODO: propagate _fflag_lp64
  _fflag_ilp64 += _fflag_preprocess
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

cython_c_args += [use_math_defines]
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
  install_tag: 'python-runtime',
)

# Ordering of subdirs: special and linalg come first, because other submodules
# have dependencies on cython_special.pxd and cython_linalg.pxd. After those,
# subdirs with the most heavy builds should come first (that parallelizes
# better)
subdir('_build_utils')
# Update BLAS/LAPACK dependencies to include wrapper libs
blas_dep = declare_dependency(
  dependencies: blas_dep,
  link_with: [g77_abi_wrappers, blas_lapack_wrapper_lib]
)
lapack_dep = declare_dependency(
  dependencies: lapack_dep,
  link_with: [g77_abi_wrappers, blas_lapack_wrapper_lib]
)

subdir('_lib')
subdir('special')
subdir('linalg')
subdir('sparse')
subdir('stats')
subdir('fft')
subdir('io')
subdir('optimize')
subdir('spatial')
subdir('cluster')
subdir('constants')
subdir('fftpack')
subdir('integrate')
subdir('differentiate')
subdir('signal')
subdir('interpolate')
subdir('ndimage')
subdir('odr')
subdir('datasets')
subdir('misc')
