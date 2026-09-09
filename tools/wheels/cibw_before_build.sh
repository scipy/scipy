set -xe

PROJECT_DIR="${1:-$PWD}"
SCIPY_SRC_DIR="${1:-$PWD}"


# Every wheel links against scipy-openblas32, except the macOS ones built against
# Accelerate - wheels.yml sets INSTALL_OPENBLAS=false for those, and leaves it unset
# everywhere else. scipy-openblas64 is not used for wheels at all.
INSTALL_OPENBLAS=${INSTALL_OPENBLAS:-true}
OPENBLAS_GRP=""
if [[ "$INSTALL_OPENBLAS" = "true" ]] ; then
    OPENBLAS_GRP="--group openblas32"
fi


# install build dependencies via uv
# log which uv this is: on Linux it comes from the manylinux image, not from setup-uv
uv --version

PYTHON_EXE="$(python -c 'import sys; print(sys.executable)')"
uv export --project "$PROJECT_DIR" --no-default-groups --group build --no-emit-project $OPENBLAS_GRP --frozen | \
    uv pip install --python "$PYTHON_EXE" --no-deps --require-hashes -r -


# Configure the pkg-config file for OpenBLAS
if [[ "$INSTALL_OPENBLAS" = "true" ]] ; then
    # The PKG_CONFIG_PATH environment variable will be pointed to this path in
    # cibuildwheel.toml and .github/workflows/wheels.yml. Note that
    # `pkgconf_path` here is only a bash variable local to this file.
    pkgconf_path=$PROJECT_DIR/.openblas
    rm -rf $pkgconf_path
    mkdir -p $pkgconf_path
    python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $pkgconf_path/scipy-openblas.pc
fi
