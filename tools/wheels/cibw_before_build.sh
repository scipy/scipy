set -xe

PROJECT_DIR="${1:-$PWD}"


if [[ $(python -c"import sys; print(sys.maxsize)") < $(python -c"import sys; print(2**33)") ]]; then
    echo "No BLAS used for 32-bit wheels"
    export INSTALL_OPENBLAS=false
elif [ -z $INSTALL_OPENBLAS ]; then
    # the macos_arm64 build might not set this variable
    export INSTALL_OPENBLAS=true
fi

if [[ $RUNNER_OS == "macOS" ]]; then
  if [ -z $INSTALL_GFORTRAN ]; then
    # the macos_arm64 build might not set this variable
    export INSTALL_GFORTRAN=true
  fi
  if [[ $INSTALL_GFORTRAN == "true" ]]; then
    source $PROJECT_DIR/tools/wheels/gfortran_macos.sh
  fi
fi

# Install OpenBLAS from scipy-openblas32|64
if [[ "$INSTALL_OPENBLAS" = "true" ]] ; then
    # By default, use scipy-openblas32
    # On 32-bit platforms and on win-arm64, use scipy-openblas32
    OPENBLAS=openblas32
    # Possible values for RUNNER_ARCH in GitHub Actions are: X86, X64, ARM, or ARM64
    if [[ $RUNNER_ARCH == "X86" || $RUNNER_ARCH == "ARM" ]] ; then
        OPENBLAS=openblas32
    elif [[ $RUNNER_ARCH == "ARM64" && $RUNNER_OS == "Windows" ]] ; then
        OPENBLAS=openblas32
    fi

    # The PKG_CONFIG_PATH environment variable will be pointed to this path in
    # cibuildwheel.toml and .github/workflows/wheels.yml. Note that
    # `pkgconf_path` here is only a bash variable local to this file.
    pkgconf_path=$PROJECT_DIR/.openblas
    echo pkgconf_path is $pkgconf_path, OPENBLAS is ${OPENBLAS}
    rm -rf $pkgconf_path
    mkdir -p $pkgconf_path
    python -m pip install -r $PROJECT_DIR/requirements/openblas.txt
    python -c "import scipy_${OPENBLAS}; print(scipy_${OPENBLAS}.get_pkg_config())" > $pkgconf_path/scipy-openblas.pc

    if [[ $RUNNER_OS == "macOS" ]]; then
      # For scipy_openblas we need the older fortran compilers that were used to
      # build it, homebrew's are too modern.

      lib_loc=$(python -c"import scipy_openblas32; print(scipy_openblas32.get_lib_dir())")
      # Use the libgfortran from gfortran rather than the one in the wheel
      # since delocate gets confused if there is more than one
      # https://github.com/scipy/scipy/issues/20852
      install_name_tool -change @loader_path/../.dylibs/libgfortran.5.dylib @rpath/libgfortran.5.dylib $lib_loc/libsci*
      install_name_tool -change @loader_path/../.dylibs/libgcc_s.1.1.dylib @rpath/libgcc_s.1.1.dylib $lib_loc/libsci*
      install_name_tool -change @loader_path/../.dylibs/libquadmath.0.dylib @rpath/libquadmath.0.dylib $lib_loc/libsci*

      codesign -s - -f $lib_loc/libsci*
    fi

    # Copy scipy-openblas DLL's to a fixed location so we can point delvewheel
    # at it in `repair_windows.sh` (needed only on Windows because of the lack
    # of RPATH support).
    if [[ $RUNNER_OS == "Windows" ]]; then
        python <<EOF
import os, scipy_${OPENBLAS}, shutil
srcdir = os.path.join(os.path.dirname(scipy_${OPENBLAS}.__file__), "lib")
shutil.copytree(srcdir, os.path.join("$pkgconf_path", "lib"))
EOF
    fi
fi

# cibuildwheel doesn't install delvewheel by default
if [[ $RUNNER_OS == "Windows" ]]; then
    python -m pip install -r $PROJECT_DIR/requirements/delvewheel_requirements.txt
    # pkgconf - carries out the role of pkg-config.
    # Alternative is pkgconfiglite that you have to install with choco
    python -m pip install -r $PROJECT_DIR/requirements/pkgconf.txt
fi
