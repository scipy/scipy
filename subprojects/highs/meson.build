project('highs', 'cpp', 'c',
  version : '1.16.0',
  meson_version: '>= 1.1.0',
  default_options : ['warning_level=1',
                     'cpp_std=c++17',
                     'wrap_mode=forcefallback',
                     'highsint64=false'])


# Add C++ compiler options
_args = [] # Extra arguments
_deps = [] # Dependencies
_incdirs = [] # All the includes

cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

# Platform detection
host_system = host_machine.system()
is_windows = host_system == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'

# Granular flags
# C++
_cpp_Wno_invalid_offsetof = cppc.get_supported_arguments('-Wno-invalid-offsetof')
_cpp_Wno_maybe_uninitialized = cppc.get_supported_arguments('-Wno-maybe-uninitialized')
_cpp_Wno_reorder = cppc.get_supported_arguments('-Wno-reorder')
_cpp_Wno_reorder_ctor = cppc.get_supported_arguments('-Wno-reorder-ctor')
_cpp_Wno_sometimes_uninitialized = cppc.get_supported_arguments('-Wno-sometimes-uninitialized')
_cpp_Wno_unused_but_set_variable = cppc.get_supported_arguments('-Wno-unused-but-set-variable')
_cpp_Wno_unused_variable = cppc.get_supported_arguments('-Wno-unused-variable')
_cpp_Wno_use_after_free = cppc.get_supported_arguments('-Wno-use-after-free')
_cpp_Wno_comment = cppc.get_supported_arguments('-Wno-comment')


# C
Wno_comment = cc.get_supported_arguments('-Wno-comment')
Wno_invalid_offsetof = cc.get_supported_arguments('-Wno-invalid-offsetof')
Wno_maybe_uninitialized = cc.get_supported_arguments('-Wno-maybe-uninitialized')
Wno_sometimes_uninitialized = cc.get_supported_arguments('-Wno-sometimes-uninitialized')
Wno_unused_label = cc.get_supported_arguments('-Wno-unused-label')
Wno_unused_but_set_variable = cc.get_supported_arguments('-Wno-unused-but-set-variable')
Wno_unused_variable = cc.get_supported_arguments('-Wno-unused-variable')
Wno_use_after_free = cc.get_supported_arguments('-Wno-use-after-free')

# Conditional arguments
_highs_cpp_suppressions = [
  _cpp_Wno_comment,
  _cpp_Wno_invalid_offsetof,
  _cpp_Wno_maybe_uninitialized,
  _cpp_Wno_reorder,
  _cpp_Wno_reorder_ctor,
  _cpp_Wno_sometimes_uninitialized,
  _cpp_Wno_unused_but_set_variable,
  _cpp_Wno_unused_variable,
  _cpp_Wno_use_after_free,
]

_highs_c_suppressions = [
  Wno_comment,
  Wno_invalid_offsetof,
  Wno_maybe_uninitialized,
  Wno_sometimes_uninitialized,
  Wno_unused_but_set_variable,
  Wno_unused_label,
  Wno_unused_variable,
  Wno_use_after_free,
]


if cppc.get_id() == 'msvc'
  add_project_arguments(
      '/wd4018',  # Disable warning: 'expression' : signed/unsigned mismatch
      '/wd4061',  # Disable warning: enumerator 'identifier' in switch of enum 'enumeration' is not explicitly handled by a case label
      '/wd4100',  # Disable warning: 'identifier' : unreferenced formal parameter
      '/wd4101',  # Disable warning: 'identifier' : unreferenced local variable
      '/wd4127',  # Disable warning: conditional expression is constant
      '/wd4189',  # Disable warning: 'identifier' : local variable is initialized but not referenced
      '/wd4244',  # Disable warning: 'conversion' conversion from 'type1' to 'type2', possible loss of data
      '/wd4245',  # Disable warning: 'conversion' conversion from 'type1' to 'type2', signed/unsigned mismatch
      '/wd4267',  # Disable warning: 'conversion' conversion from 'size_t' to 'type', possible loss of data
      '/wd4324',  # Disable warning: 'structure' structure was padded due to alignment specifier
      '/wd4365',  # Disable warning: 'action' : conversion from 'type_1' to 'type_2', signed/unsigned mismatch
      '/wd4389',  # Disable warning: 'modifier' : signed/unsigned mismatch
      '/wd4456',  # Disable warning: declaration of 'identifier' hides previous local declaration
      '/wd4457',  # Disable warning: declaration of 'identifier' hides function parameter
      '/wd4458',  # Disable warning: declaration of 'identifier' hides class member
      '/wd4459',  # Disable warning: declaration of 'identifier' hides global declaration
      '/wd4514',  # Disable warning: 'function' : unreferenced inline function has been removed
      '/wd4701',  # Disable warning: potentially uninitialized local variable 'name' used
      '/wd4820',  # Disable warning: 'bytes' bytes padding added after construct 'member_name'
      language: 'cpp',
  )
  _args += '-D_CRT_SECURE_NO_WARNINGS'
endif

cpu_family = host_machine.cpu_family()

if cpu_family in ['x86_64', 'i686'] and not is_windows
    add_project_arguments(cppc.get_supported_arguments('-mpopcnt'), language: 'cpp')
endif
if cpu_family in ['ppc64', 'powerpc64'] and not meson.is_cross_build()
    add_project_arguments(cppc.get_supported_arguments('-mpopcntd'), language: 'cpp')
endif

if is_mingw
  # For mingw-w64, don't use LTO
  add_project_arguments('-fno-use-linker-plugin', language: ['c', 'cpp'])
endif

# --------------------- Dependencies
threads_dep = dependency('threads', required: false)
_deps += threads_dep

# Determine whether it is necessary to link libatomic. This could be the case
# e.g. on 32-bit platforms when atomic operations are used on 64-bit types.
# The check is copied from SciPy which in turn came from
# Mesa <https://www.mesa3d.org/>.
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
      # From SciPy
      # We're not sure that with `-latomic` things will work for all compilers,
      # so verify and only keep libatomic as a dependency if this works.
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

_deps += atomic_dep

# Optional
zlib_dep = dependency('zlib',
                      required: get_option('use_zlib'))
if zlib_dep.found()
  _deps += zlib_dep
  _incdirs += include_directories(['extern/zstr'])
endif

# Include Sources
_incdirs += include_directories([
  'extern/',
  'extern/catch',
  'extern/pdqsort/',
])

# Optional arguments
if get_option('fast_build')
  add_project_arguments(cc.get_supported_arguments('-fno-omit-frame-pointer'), language : ['c', 'cpp'])
endif

if get_option('debug_sol')
  _args += ['-DHIGHS_DEBUGSOL']
endif

if cppc.get_id() == 'msvc'
  # Don't depend on VCRUNTIME140_1.dll
  # https://cibuildwheel.readthedocs.io/en/stable/faq/#windows-importerror-dll-load-failed-the-specific-module-could-not-be-found
  add_project_arguments('/d2FH4-', language : ['cpp', 'c'])
  if get_option('with_stdcall')
    add_project_arguments(cppc.get_supported_arguments('/Gz'), language : ['cpp', 'c'])
  endif
endif

# --------------------- Library

subdir('highs') # defines highslib

# --------------------- Tests

if get_option('with_tests')
  subdir('check')
  _deps += dependency('catch2', required: true)
endif

# --------------------- Bindings
# now in highs/
