set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")
echo $PLATFORM

# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_osx.txt >> $PROJECT_DIR/LICENSE.txt

#########################################################################################
# Install GFortran + OpenBLAS

if [[ $PLATFORM == "macosx-x86_64" ]]; then
  # Openblas
  basedir=$(python tools/openblas_support.py)

  # copy over the OpenBLAS library stuff first
  cp -r $basedir/lib/* /usr/local/lib
  cp $basedir/include/* /usr/local/include

  #GFORTRAN=$(type -p gfortran-9)
  #sudo ln -s $GFORTRAN /usr/local/bin/gfortran
  # same version of gfortran as the openblas-libs and scipy-wheel builds
  curl -L https://github.com/isuruf/gcc/releases/download/gcc-11.3.0-2/gfortran-darwin-x86_64-native.tar.gz -o gfortran.tar.gz

  GFORTRAN_SHA256=$(shasum -a 256 gfortran.tar.gz)
  KNOWN_SHA256="981367dd0ad4335613e91bbee453d60b6669f5d7e976d18c7bdb7f1966f26ae4  gfortran.tar.gz"
  if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
      echo sha256 mismatch
      exit 1
  fi

  sudo mkdir -p /opt/
  # places gfortran in /opt/gfortran-darwin-x86_64-native. There's then
  # bin, lib, include, libexec underneath that.
  sudo tar -xv -C /opt -f gfortran.tar.gz

  # Link these into /usr/local so that there's no need to add rpath or -L
  for f in libgfortran.dylib libgfortran.5.dylib libgcc_s.1.dylib libgcc_s.1.1.dylib libquadmath.dylib libquadmath.0.dylib; do
    ln -sf /opt/gfortran-darwin-x86_64-native/lib/$f /usr/local/lib/$f
  done
  ln -sf /opt/gfortran-darwin-x86_64-native/bin/gfortran /usr/local/bin/gfortran

  # Set SDKROOT env variable if not set
  # This step is required whenever the gfortran compilers sourced from
  # conda-forge (built by isuru fernando) are used outside of a conda-forge
  # environment (so it mirrors what is done in the conda-forge compiler
  # activation scripts)
  export SDKROOT=${SDKROOT:-$(xcrun --show-sdk-path)}
fi

if [[ $PLATFORM == "macosx-arm64" ]]; then
  # OpenBLAS
  # need a version of OpenBLAS that is suited for gcc >= 11
  basedir=$(python tools/openblas_support.py)

  # use /opt/arm64-builds as a prefix, because that's what the multibuild
  # OpenBLAS pkgconfig files state
  sudo mkdir -p /opt/arm64-builds/lib
  sudo mkdir -p /opt/arm64-builds/include
  sudo cp -r $basedir/lib/* /opt/arm64-builds/lib
  sudo cp $basedir/include/* /opt/arm64-builds/include

  # we want to force a dynamic linking
  sudo rm /opt/arm64-builds/lib/*.a

  curl -L https://github.com/fxcoudert/gfortran-for-macOS/releases/download/12.1-monterey/gfortran-ARM-12.1-Monterey.dmg -o gfortran.dmg
  GFORTRAN_SHA256=$(shasum -a 256 gfortran.dmg)
  KNOWN_SHA256="e2e32f491303a00092921baebac7ffb7ae98de4ca82ebbe9e6a866dd8501acdf  gfortran.dmg"

  if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
      echo sha256 mismatch
      exit 1
  fi

  hdiutil attach -mountpoint /Volumes/gfortran gfortran.dmg
  sudo installer -pkg /Volumes/gfortran/gfortran.pkg -target /
  type -p gfortran
fi
