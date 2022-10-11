project(
  'SciPy',
  'c', 'cpp',
  # Note that the git commit hash cannot be added dynamically here (it is added
  # in the dynamically generated and installed `scipy/version.py` though - see
  # tools/version_utils.py
  version: '1.10.0.dev0',
  license: 'BSD-3',
  meson_version: '>= 0.63.0',
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

# Check compiler is recent enough (see "Toolchain Roadmap" for details)
if cc.get_id() == 'gcc'
  if not cc.version().version_compare('>=8.0')
    error('SciPy requires GCC >= 8.0')
  endif
elif cc.get_id() == 'msvc'
  if not cc.version().version_compare('>=19.20')
    error('SciPy requires at least vc142 (default with Visual Studio 2019) ' + \
          'when building with MSVC')
  endif
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
if host_machine.system() == 'os400'
  add_global_arguments('-pthread', language : 'cpp')
  add_global_arguments('-D__STDC_FORMAT_MACROS', language : 'cpp')
  add_project_link_arguments('-Wl,-bnotextro', language : 'c')
  add_project_link_arguments('-Wl,-bnotextro', language : 'cpp')
  add_project_link_arguments('-Wl,-bnotextro', language : 'fortran')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
if ff.has_argument('-Wno-conversion')
  add_project_arguments('-Wno-conversion', language: 'fortran')
endif

is_windows = host_machine.system() == 'windows'

# Intel compilers default to fast-math, so disable it if we detect Intel
# compilers. A word of warning: this may not work with the conda-forge
# compilers, because those have the annoying habit of including lots of flags
# that are gcc-specific in CFLAGS/CXXFLAGS/FFLAGS, which throws off the
# detection logic below. You have to remove the wrong flags (only `-isystem`
# is actually needed, everything else shouldn't be there).
_intel_cflags = []
_intel_fflags = []
if cc.get_id() == 'intel'
  _intel_cflags += cc.get_supported_arguments('-fp-model=strict')
elif cc.get_id() == 'intel-cl'
  _intel_cflags += cc.get_supported_arguments('/fp:strict')
endif
if ff.get_id() == 'intel'
  _intel_fflags = ff.get_supported_arguments('-fp-model=strict')
  minus0_arg = ['-assume', 'minus0']
  if ff.has_multi_arguments(minus0_arg)
    _intel_fflags += minus0_arg
  endif
elif ff.get_id() == 'intel-cl'
  # Intel Fortran on Windows does things differently, so deal with that
  # (also specify dynamic linking and the right name mangling)
  _intel_fflags = ff.get_supported_arguments(
    '/fp:strict', '/MD', '/names:lowercase', '/assume:underscore',
    '/assume:minus0'
  )
endif
add_project_arguments(_intel_cflags, language: ['c', 'cpp'])
add_project_arguments(_intel_fflags, language: 'fortran')

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
py3_dep = py3.dependency()

subdir('scipy')
