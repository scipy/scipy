# shellcheck shell=sh

# This function takes no arguments
# It tries to determine the name of this file in a programatic way.
_get_sourced_filename() {
    # shellcheck disable=SC3054,SC2296 # non-POSIX array access and bad '(' are guarded
    if [ -n "${BASH_SOURCE+x}" ] && [ -n "${BASH_SOURCE[0]}" ]; then
        # shellcheck disable=SC3054 # non-POSIX array access is guarded
        basename "${BASH_SOURCE[0]}"
    elif [ -n "$ZSH_NAME" ] && [ -n "${(%):-%x}" ]; then
        # in zsh use prompt-style expansion to introspect the same information
        # see http://stackoverflow.com/questions/9901210/bash-source0-equivalent-in-zsh
        # shellcheck disable=SC2296  # bad '(' is guarded
        basename "${(%):-%x}"
    else
        echo "UNKNOWN FILE"
    fi
}

# The arguments to this are:
# 1. activation nature {activate|deactivate}
# 2. toolchain nature {build|host|ccc}
# 3. machine (should match -dumpmachine)
# 4. prefix (including any final -)
# 5+ program (or environment var comma value)
# The format for 5+ is name{,,value}. If value is specified
#  then name taken to be an environment variable, otherwise
#  it is taken to be a program. In this case, which is used
#  to find the full filename during activation. The original
#  value is stored in environment variable CONDA_BACKUP_NAME
#  For deactivation, the distinction is irrelevant as in all
#  cases NAME simply gets reset to CONDA_BACKUP_NAME.  It is
#  a fatal error if a program is identified but not present.
_tc_activation() {
  local act_nature="$1"; shift
  local tc_prefix="$1"; shift
  local thing
  local newval
  local from
  local to
  local pass

  if [ "${act_nature}" = "activate" ]; then
    from=""
    to="CONDA_BACKUP_"
  else
    from="CONDA_BACKUP_"
    to=""
  fi

  for pass in check apply; do
    for thing in "$@"; do
      case "${thing}" in
        *,*)
          newval=$(echo "${thing}" | sed "s,^[^\,]*\,\(.*\),\1,")
          thing=$(echo "${thing}" | sed "s,^\([^\,]*\)\,.*,\1,")
          ;;
        *)
          newval="${CONDA_PREFIX}/bin/${tc_prefix}${thing}"
          thing=$(echo "${thing}" | tr 'a-z+-' 'A-ZX_')
          if [ ! -x "${newval}" ] && [ "${pass}" = "check" ]; then
            echo "ERROR: This cross-compiler package contains no program ${newval}"
            return 1
          fi
          ;;
      esac
      if [ "${pass}" = "apply" ]; then
        eval oldval="\$${from}$thing"
        if [ -n "${oldval}" ]; then
          eval export "${to}'${thing}'=\"${oldval}\""
        else
          eval unset '${to}${thing}'
        fi
        if [ -n "${newval}" ]; then
          eval export "'${from}${thing}=${newval}'"
        else
          eval unset '${from}${thing}'
        fi
      fi
    done
  done
  return 0
}

# The compiler adds $PREFIX/lib to rpath, so it's better to add -L and -isystem  as well.
if [ "${CONDA_BUILD:-0}" = "1" ]; then
  CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${PREFIX}/include -fdebug-prefix-map=${SRC_DIR}=/usr/local/src/conda/${PKG_NAME}-${PKG_VERSION} -fdebug-prefix-map=${PREFIX}=/usr/local/src/conda-prefix"
  DEBUG_CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem ${PREFIX}/include -fdebug-prefix-map=${SRC_DIR}=/usr/local/src/conda/${PKG_NAME}-${PKG_VERSION} -fdebug-prefix-map=${PREFIX}=/usr/local/src/conda-prefix"
  LDFLAGS_USED="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,${PREFIX}/lib -Wl,-rpath-link,${PREFIX}/lib -L${PREFIX}/lib"
  CPPFLAGS_USED="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem ${PREFIX}/include"
  DEBUG_CPPFLAGS_USED="-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem ${PREFIX}/include"
  CMAKE_PREFIX_PATH_USED="${PREFIX}:${CONDA_PREFIX}/x86_64-conda-linux-gnu/sysroot/usr"
else
  CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem ${CONDA_PREFIX}/include"
  DEBUG_CFLAGS_USED="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -isystem ${CONDA_PREFIX}/include"
  CPPFLAGS_USED="-DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem ${CONDA_PREFIX}/include"
  DEBUG_CPPFLAGS_USED="-D_DEBUG -D_FORTIFY_SOURCE=2 -Og -isystem ${CONDA_PREFIX}/include"
  LDFLAGS_USED="-Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,${CONDA_PREFIX}/lib -Wl,-rpath-link,${CONDA_PREFIX}/lib -L${CONDA_PREFIX}/lib"
  CMAKE_PREFIX_PATH_USED="${CONDA_PREFIX}:${CONDA_PREFIX}/x86_64-conda-linux-gnu/sysroot/usr"
fi

if [ "${CONDA_BUILD:-0}" = "1" ]; then
  if [ -f /tmp/old-env-$$.txt ]; then
    rm -f /tmp/old-env-$$.txt || true
  fi
  env > /tmp/old-env-$$.txt
fi

_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED=${_CONDA_PYTHON_SYSCONFIGDATA_NAME:-_sysconfigdata_x86_64_conda_cos6_linux_gnu}
if [ -n "${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}" ] && [ -n "${SYS_SYSROOT}" ]; then
  if find "$(dirname "$(dirname "${SYS_PYTHON}")")/lib/"python* -type f -name "${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}.py" -exec false {} +; then
    echo ""
    echo "WARNING: The Python interpreter at the following prefix:"
    echo "         $(dirname "$(dirname "${SYS_PYTHON}")")"
    echo "         .. is not able to handle sysconfigdata-based compilation for the host:"
    echo "         $( printf %s "${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}" | sed s/_sysconfigdata_//g )"
    echo ""
    echo "         We are not preventing things from continuing here, but *this* Python will not"
    echo "         be able to compile software for this host, and, depending on whether it has"
    echo "         been patched to ignore missing _CONDA_PYTHON_SYSCONFIGDATA_NAME or not, may cause"
    echo "         an exception."
    echo ""
    echo "         This can happen for one of three reasons:"
    echo ""
    echo "         1. It is out of date: Please run 'conda update python' in that environment"
    echo ""
    echo "         2. You are bootstrapping a sysconfigdata-based cross-capable Python and can ignore this"
    echo "            (but please remember to copy the generated sysconfigdata back to the Python recipe's"
    echo "             sysconfigdate folder and then rebuild it for all the systems you want to be able"
    echo "             to use as a build machine for this host)."
    echo ""
    echo "         3. You are attempting your own bespoke cross-compilation host that is not supported. Have"
    echo "            you provided your own value in the _CONDA_PYTHON_SYSCONFIGDATA_NAME environment variable but"
    echo "            misspelt it and/or failed to add the neccessary ${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}.py"
    echo "            file to the Python interpreter's standard library?"
    echo ""
  fi
fi

_CMAKE_ARGS="-DCMAKE_AR=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-ar -DCMAKE_CXX_COMPILER_AR=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc-ar -DCMAKE_C_COMPILER_AR=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc-ar"
_CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_RANLIB=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-ranlib -DCMAKE_CXX_COMPILER_RANLIB=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc-ranlib -DCMAKE_C_COMPILER_RANLIB=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gcc-ranlib"
_CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_LINKER=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-ld -DCMAKE_STRIP=${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-strip"
_CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_BUILD_TYPE=Release"

_MESON_ARGS="--buildtype release"

if [ "${CONDA_BUILD:-0}" = "1" ]; then
  _CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY"
  _CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_FIND_ROOT_PATH=$PREFIX;${BUILD_PREFIX}/x86_64-conda-linux-gnu/sysroot"
  _CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_INSTALL_LIBDIR=lib"
  _CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_PROGRAM_PATH=${BUILD_PREFIX}/bin;$PREFIX/bin"
  _MESON_ARGS="${_MESON_ARGS} --prefix=$PREFIX -Dlibdir=lib"
fi

# shellcheck disable=SC2050 # templating will fix this error
if [ "" = "1" ]; then
  _CMAKE_ARGS="${_CMAKE_ARGS} -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_SYSTEM_PROCESSOR=x86_64"
  _MESON_ARGS="${_MESON_ARGS} --cross-file ${CONDA_PREFIX}/meson_cross_file.txt"
  echo "[host_machine]" > "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "system = 'linux'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "cpu = 'x86_64'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "cpu_family = 'x86_64'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "endian = 'little'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  # specify path to correct binaries from build (not host) environment,
  # which meson will not auto-discover (out of caution) if not told explicitly.
  echo "[binaries]" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "cmake = '${CONDA_PREFIX}/bin/cmake'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
  echo "pkgconfig = '${CONDA_PREFIX}/bin/pkg-config'" >> "${CONDA_PREFIX}/meson_cross_file.txt"
fi

_tc_activation \
  activate x86_64-conda-linux-gnu- "HOST,x86_64-conda-linux-gnu" "BUILD,x86_64-conda-linux-gnu" \
  "CONDA_TOOLCHAIN_HOST,x86_64-conda-linux-gnu" \
  "CONDA_TOOLCHAIN_BUILD,x86_64-conda-linux-gnu" \
  cc cpp gcc gcc-ar gcc-nm gcc-ranlib \
  "CPPFLAGS,${CPPFLAGS_USED}${CPPFLAGS:+ }${CPPFLAGS:-}" \
  "CFLAGS,${CFLAGS_USED}${CFLAGS:+ }${CFLAGS:-}" \
  "LDFLAGS,${LDFLAGS_USED}${LDFLAGS:+ }${LDFLAGS:-}" \
  "DEBUG_CPPFLAGS,${DEBUG_CPPFLAGS_USED}${DEBUG_CPPFLAGS:+ }${DEBUG_CPPFLAGS:-}" \
  "DEBUG_CFLAGS,${DEBUG_CFLAGS_USED}${DEBUG_CFLAGS:+ }${DEBUG_CFLAGS:-}" \
  "CMAKE_PREFIX_PATH,${CMAKE_PREFIX_PATH_USED}" \
  "_CONDA_PYTHON_SYSCONFIGDATA_NAME,${_CONDA_PYTHON_SYSCONFIGDATA_NAME_USED}" \
  "CONDA_BUILD_SYSROOT,${CONDA_PREFIX}/x86_64-conda-linux-gnu/sysroot" \
  "CONDA_BUILD_CROSS_COMPILATION," \
  "CC_FOR_BUILD,${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-cc" \
  "build_alias,x86_64-conda-linux-gnu" \
  "host_alias,x86_64-conda-linux-gnu" \
  "MESON_ARGS,${_MESON_ARGS}" \
  "CMAKE_ARGS,${_CMAKE_ARGS}"

unset _CMAKE_ARGS

# shellcheck disable=SC2050 # templating will fix this error
if [ "" = "1" ]; then
_tc_activation \
   activate x86_64-conda-linux-gnu- \
   "QEMU_LD_PREFIX,${QEMU_LD_PREFIX:-${CONDA_BUILD_SYSROOT}}"
fi

if [ $? -ne 0 ]; then
  echo "ERROR: $(_get_sourced_filename) failed, see above for details"
else
  if [ "${CONDA_BUILD:-0}" = "1" ]; then
    if [ -f /tmp/new-env-$$.txt ]; then
      rm -f /tmp/new-env-$$.txt || true
    fi
    env > /tmp/new-env-$$.txt

    echo "INFO: $(_get_sourced_filename) made the following environmental changes:"
    diff -U 0 -rN /tmp/old-env-$$.txt /tmp/new-env-$$.txt | tail -n +4 | grep "^-.*\|^+.*" | grep -v "CONDA_BACKUP_" | sort
    rm -f /tmp/old-env-$$.txt /tmp/new-env-$$.txt || true
  fi

  # fix prompt for zsh
  if [ -n "${ZSH_NAME:-}" ]; then
    autoload -Uz add-zsh-hook

    _conda_clang_precmd() {
      # shellcheck disable=SC2034 # guarded for zsh
      HOST="${CONDA_BACKUP_HOST}"
    }
    add-zsh-hook -Uz precmd _conda_clang_precmd

    _conda_clang_preexec() {
      # shellcheck disable=SC2034 # guarded for zsh
      HOST="${CONDA_TOOLCHAIN_HOST}"
    }
    add-zsh-hook -Uz preexec _conda_clang_preexec
  fi
fi
