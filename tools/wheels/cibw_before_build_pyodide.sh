#!/bin/bash
# Provision BLAS/LAPACK for SciPy's WebAssembly/Pyodide build.
#
# Unlike the native builds (which install scipy-openblas via
# tools/wheels/cibw_before_build.sh), the WASM build does not
# use OpenBLAS. Instead we build, with Emscripten:
#
#   1. BLIS               - The BLAS backend, from https://github.com/pyodide/blis,
#                           which is a fork of https://github.com/flame/blis
#   2. semicolon-lapack   - A pure C reimplementation of LAPACK, available at
#                           https://github.com/ilayn/semicolon-lapack
#
#
# For more info, please see:
# - https://github.com/pyodide/pyodide-recipes/issues/604
# - https://github.com/pyodide/pyodide-recipes/pull/619

set -xeo pipefail

PROJECT_DIR="${1:-$PWD}"
PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"

# cibuildwheel doesn't have a default wheel repair command for Pyodide
python -m pip install "auditwheel-emscripten==0.2.5"

# for semicolon-lapack to build
python -m pip install "meson>=1.5.0" "ninja>=1.8.2"

# BLIS
BLIS_REPO="https://github.com/pyodide/blis"
BLIS_REF="pyodide-2.1"

# semicolon-lapack
SEMILAPACK_VERSION="0.01.3-pre"
SEMILAPACK_URL="https://github.com/ilayn/semicolon-lapack/archive/refs/tags/v${SEMILAPACK_VERSION}.tar.gz"

# To be consumed via PKG_CONFIG_PATH
BLAS_PREFIX="${PROJECT_DIR}/.blis"

BUILD_ROOT="$(mktemp -d)"
trap 'rm -rf "${BUILD_ROOT}"' EXIT

rm -rf "${BLAS_PREFIX}"
mkdir -p "${BLAS_PREFIX}"

# cibuildwheel runs before-build with PYODIDE_ROOT set, but not the individual
# build variables (those are injected by pyodide-build during the wheel build
# itself). Recover them here through `pyodide config get`, which reads the
# xbuildenv configuration, so we build BLIS/semicolon-lapack with exactly the
# same Emscripten cross file and linker flags SciPy will be built with.
MESON_CROSS_FILE="${MESON_CROSS_FILE:-$(pyodide config get meson_cross_file)}"
SIDE_MODULE_LDFLAGS="${SIDE_MODULE_LDFLAGS:-$(pyodide config get ldflags)}"

# N.B. BLIS and semicolon-lapack install their .pc files under different
# prefixes (lib/pkgconfig and share/pkgconfig), so we have to search both.
export PKG_CONFIG_PATH="${BLAS_PREFIX}/lib/pkgconfig:${BLAS_PREFIX}/share/pkgconfig:${PKG_CONFIG_PATH:-}"

# ---------------------------------------------------------------------------
# 1. BLIS
# ---------------------------------------------------------------------------

git clone --depth 1 --branch "${BLIS_REF}" "${BLIS_REPO}" "${BUILD_ROOT}/blis"
pushd "${BUILD_ROOT}/blis"
    CC=emcc AR=emar RANLIB=emranlib \
        ./configure \
            --disable-threading \
            --disable-shared \
            --enable-cblas \
            --prefix="${BLAS_PREFIX}" \
            wasm32
    emmake make -j "${PYODIDE_JOBS:-3}"
    emmake make install

    emcc "${BLAS_PREFIX}/lib/libblis.a" \
        ${SIDE_MODULE_LDFLAGS} \
        -o "${BLAS_PREFIX}/lib/libblis.so"
popd

# ---------------------------------------------------------------------------
# 2. semicolon-lapack
# ---------------------------------------------------------------------------

curl -L "${SEMILAPACK_URL}" -o "${BUILD_ROOT}/semilapack.tar.gz"
mkdir -p "${BUILD_ROOT}/semilapack"
tar -xzf "${BUILD_ROOT}/semilapack.tar.gz" -C "${BUILD_ROOT}/semilapack" --strip-components=1
pushd "${BUILD_ROOT}/semilapack"
    CC=emcc CXX=em++ AR=emar \
        meson setup builddir \
            --cross-file "${MESON_CROSS_FILE}" \
            -Dblas=blis \
            -Dfabi_shim=true \
            -Dtests=false \
            --default-library=static \
            --prefix="${BLAS_PREFIX}"
    ninja -C builddir
    ninja -C builddir install

    # As with BLIS, link a shared side module, resolving CBLAS symbols against
    # libblis.so.
    emcc "${BLAS_PREFIX}/lib/libsemilapack.a" \
        "${BLAS_PREFIX}/lib/libsemilapack_fortran.a" \
        "${BLAS_PREFIX}/lib/libblis.so" \
        ${SIDE_MODULE_LDFLAGS} \
        -o "${BLAS_PREFIX}/lib/libsemilapack.so"
popd

echo "WASM BLAS/LAPACK compiled in ${BLAS_PREFIX}:"
ls -l "${BLAS_PREFIX}/lib" || true
ls -l "${BLAS_PREFIX}/lib/pkgconfig" || true
