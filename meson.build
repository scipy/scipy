project(
  'scipy',
  'c', 'cpp', 'cython',
  version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
  license: 'BSD-3',
  meson_version: '>= 1.5.0',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c17',
    'cpp_std=c++17',
    'blas=openblas',
    'lapack=openblas'
  ],
)

meson.add_dist_script('tools/trim_sdist_content.py')

py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

min_numpy_version = '1.25.2'  # keep in sync with pyproject.toml
min_python_version = '3.11'   # keep in sync with pyproject.toml

python_version = py3.language_version()
if python_version.version_compare(f'<@min_python_version@')
  error(f'Minimum supported Python version is @min_python_version@, found @python_version@')
endif

# Emit a warning for 32-bit Python installs on Windows; users are getting
# unexpected from-source builds there because we no longer provide wheels.
is_windows = host_machine.system() == 'windows'
if is_windows and py3.has_variable('EXT_SUFFIX')
  ext_suffix = py3.get_variable('EXT_SUFFIX')
  if ext_suffix.contains('win32')
    warning('You are building from source on a 32-bit Python install. SciPy does not provide 32-bit wheels; install 64-bit Python if you are having issues!')
  endif
endif

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
cy = meson.get_compiler('cython')
# generator() doesn't accept compilers, only found programs - cast it.
cython = find_program(cy.cmd_array()[0])

# Check compiler is recent enough (see "Toolchain Roadmap" for details)
if cc.get_id() == 'gcc'
  if not cc.version().version_compare('>=9.1')
    error('SciPy requires GCC >= 9.1')
  endif
elif cc.get_id() == 'clang' or cc.get_id() == 'clang-cl'
  if not cc.version().version_compare('>=15.0')
    error('SciPy requires clang >= 15.0')
  endif
elif cc.get_id() == 'msvc'
  if not cc.version().version_compare('>=19.20')
    error('SciPy requires at least vc142 (default with Visual Studio 2019) ' + \
          'when building with MSVC')
  endif
endif
if not cy.version().version_compare('>=3.0.8')
  error('SciPy requires Cython >= 3.0.8')
endif

_global_c_args = cc.get_supported_arguments(
  '-Wno-unused-but-set-variable',
  '-Wno-unused-function',
  '-Wno-conversion',
  '-Wno-misleading-indentation',
)
add_project_arguments(_global_c_args, language : 'c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

if host_machine.system() == 'os400'
  # IBM i system, needed to avoid build errors - see gh-17193
  add_project_arguments('-D__STDC_FORMAT_MACROS', language : 'cpp')
  add_project_link_arguments('-Wl,-bnotextro', language : ['c', 'cpp', 'fortran'])
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.get_id() == 'gcc'
  # -std=legacy is not supported by all Fortran compilers, but very useful with
  # gfortran since it avoids a ton of warnings that we don't care about.
  # Needs fixing in Meson, see https://github.com/mesonbuild/meson/issues/11633.
  add_project_arguments('-std=legacy', language: 'fortran')
endif

if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

if host_machine.system() == 'darwin'
  if cc.has_link_argument('-Wl,-dead_strip')
    # Allow linker to strip unused symbols
    add_project_link_arguments('-Wl,-dead_strip', language : ['c', 'cpp', 'fortran'])
  endif
endif

# Intel compilers default to fast-math, so disable it if we detect Intel
# compilers. A word of warning: this may not work with the conda-forge
# compilers, because those have the annoying habit of including lots of flags
# that are gcc-specific in CFLAGS/CXXFLAGS/FFLAGS, which throws off the
# detection logic below. You have to remove the wrong flags (only `-isystem`
# is actually needed, everything else shouldn't be there).
_intel_cflags = []
_intel_fflags = []
if cc.get_id() in ['intel', 'intel-llvm']
  _intel_cflags += cc.get_supported_arguments('-fp-model=strict')
elif cc.get_id() in ['intel-cl', 'intel-llvm-cl']
  _intel_cflags += cc.get_supported_arguments('/fp:strict')
endif
if ff.get_id() in ['intel', 'intel-llvm']
  _intel_fflags = ff.get_supported_arguments('-fp-model=strict')
  minus0_arg = ['-assume', 'minus0']
  if ff.has_multi_arguments(minus0_arg)
    _intel_fflags += minus0_arg
  endif
elif ff.get_id() in ['intel-cl', 'intel-llvm-cl']
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _intel_fflags = ff.get_supported_arguments(
    '/fp:strict', '/MD', '/names:lowercase', '/assume:underscore',
    '/assume:minus0'
  )
endif
add_global_arguments(_intel_cflags, language: ['c', 'cpp'])
add_global_arguments(_intel_fflags, language: 'fortran')

# Hide symbols when building on Linux with GCC. For Python extension modules,
# we only need `PyInit_*` to be public, anything else may cause problems. So we
# use a linker script to avoid exporting those symbols (this is in addition to
# Meson using `-fvisibility=hidden` for C and `-fvisibility-inlines-hidden` for
# C++ code. See gh-15996 for details.
_linker_script = meson.project_source_root() / 'scipy/_build_utils/link-version-pyinit.map'
version_link_args = ['-Wl,--version-script=' + _linker_script]
# Note that FreeBSD only accepts version scripts when -shared is passed,
# hence we need to pass that to `cc.links` explicitly (flag is already
# present for `extension_module` invocations).
if not cc.links('', name: '-Wl,--version-script', args: ['-shared', version_link_args])
  version_link_args = []
endif

generate_f2pymod = find_program('tools/generate_f2pymod.py')
tempita = find_program('scipy/_build_utils/tempita.py')

use_pythran = get_option('use-pythran')
if use_pythran
  pythran = find_program('pythran', native: true, version: '>=0.14.0')
  # xsimd is unvendored from pythran by conda-forge, and due to a compiler
  # activation bug the default <prefix>/include/ may not be visible (see
  # gh-15698). Hence look for xsimd explicitly.
  xsimd_dep = dependency('xsimd', required: false)
endif

fs = import('fs')
if not fs.exists('subprojects/xsf/README.md')
  error('Missing the `xsf` submodule! Run `git submodule update --init` to fix this.')
endif

xsf = subproject('xsf')
xsf_dep = xsf.get_variable('xsf_dep')

use_system_libraries = get_option('use-system-libraries')
all_system_libraries = false
auto_system_libraries = false

if use_system_libraries.contains('none')
  use_system_libraries = ['none']
elif use_system_libraries.contains('all')
  all_system_libraries = true
elif use_system_libraries.contains('auto')
  auto_system_libraries = true 
endif

if all_system_libraries or use_system_libraries.contains('boost.math')
  boost_math_dep = dependency('boost', version : '1.88.0')
elif auto_system_libraries
  boost_math_dep = dependency(
    'boost', version : '1.88.0',
    fallback : ['boost_math', 'boost_math_dep']
  )
else
  boost_math = subproject('boost_math', version : '1.88.0')
  boost_math_dep = boost_math.get_variable('boost_math_dep')
endif

if all_system_libraries or use_system_libraries.contains('qhull')
  qhull_r_dep = dependency('qhull_r', version : '8.0.2')
elif auto_system_libraries
  qhull_r_dep = dependency(
    'qhull_r', version : '8.0.2',
    fallback : ['qhull_r', 'qhull_r_dep']
  )
else
  qhull_r = subproject('qhull_r', version : '8.0.2')
  qhull_r_dep = qhull_r.get_variable('qhull_r_dep')
endif

subdir('scipy')
