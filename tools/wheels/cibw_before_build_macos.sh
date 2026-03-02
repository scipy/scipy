set -xe

PROJECT_DIR="$1"
PLATFORM=$(uname -m)
echo $PLATFORM


#########################################################################################
# Install GFortran + OpenBLAS

if [[ $PLATFORM == "x86_64" ]]; then
  #GFORTRAN=$(type -p gfortran-9)
  #sudo ln -s $GFORTRAN /usr/local/bin/gfortran
  # same version of gfortran as the openblas-libs
  # https://github.com/MacPython/gfortran-install.git
  curl -L https://github.com/isuruf/gcc/releases/download/gcc-15.2.0/gfortran-darwin-x86_64-native.tar.gz -o gfortran.tar.gz
  KNOWN_SHA256="fb03c1f37bf0258ada6e3e41698e3ad416fff4dad448fd746e01d8ccf1efdc0f  gfortran.tar.gz"
  GFORTRAN_LOC=/opt/gfortran-darwin-x86_64-native
fi
if [[ $PLATFORM == "arm64" ]]; then
  curl -L https://github.com/isuruf/gcc/releases/download/gcc-15.2.0/gfortran-darwin-arm64-native.tar.gz -o gfortran.tar.gz
  KNOWN_SHA256="999a91eef894d32f99e3b641520bef9f475055067f301f0f1947b8b716b5922a  gfortran.tar.gz"
  GFORTRAN_LOC=/opt/gfortran-darwin-arm64-native
fi

  GFORTRAN_SHA256=$(shasum -a 256 gfortran.tar.gz)
  if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
      echo sha256 mismatch
      exit 1
  fi

  sudo mkdir -p /opt/
  # the following directories may not exist on a default image
  sudo mkdir -p /usr/local/bin
  sudo mkdir -p /usr/local/lib
  # places gfortran in /opt/gfortran-darwin-x86_64-native. There's then
  # bin, lib, include, libexec underneath that.
  sudo tar -xv -C /opt -f gfortran.tar.gz

  # Link these into /usr/local so that there's no need to add rpath or -L
  for f in libgfortran.dylib libgfortran.5.dylib libgcc_s.1.dylib libgcc_s.1.1.dylib libquadmath.dylib libquadmath.0.dylib; do
    sudo ln -sf $GFORTRAN_LOC/lib/$f /usr/local/lib/$f
  done
  sudo ln -sf $GFORTRAN_LOC/bin/gfortran /usr/local/bin/gfortran

  # Set SDKROOT env variable if not set
  # This step is required whenever the gfortran compilers sourced from
  # conda-forge (built by isuru fernando) are used outside of a conda-forge
  # environment (so it mirrors what is done in the conda-forge compiler
  # activation scripts)
  export SDKROOT=${SDKROOT:-$(xcrun --show-sdk-path)}

# Install OpenBLAS
python -m pip install -r requirements/openblas.txt
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PROJECT_DIR/scipy-openblas.pc

lib_loc=$(python -c"import scipy_openblas32; print(scipy_openblas32.get_lib_dir())")
# Use the libgfortran from gfortran rather than the one in the wheel
# since delocate gets confused if there is more than one
# https://github.com/scipy/scipy/issues/20852
install_name_tool -change @loader_path/../.dylibs/libgfortran.5.dylib @rpath/libgfortran.5.dylib $lib_loc/libsci*
install_name_tool -change @loader_path/../.dylibs/libgcc_s.1.1.dylib @rpath/libgcc_s.1.1.dylib $lib_loc/libsci*
install_name_tool -change @loader_path/../.dylibs/libquadmath.0.dylib @rpath/libquadmath.0.dylib $lib_loc/libsci*

codesign -s - -f $lib_loc/libsci*
