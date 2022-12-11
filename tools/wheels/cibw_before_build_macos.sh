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
  curl -L https://github.com/MacPython/gfortran-install/raw/master/archives/gfortran-4.9.0-Mavericks.dmg -o gfortran.dmg
  GFORTRAN_SHA256=$(shasum -a 256 gfortran.dmg)
  KNOWN_SHA256="d2d5ca5ba8332d63bbe23a07201c4a0a5d7e09ee56f0298a96775f928c3c4b30  gfortran.dmg"
  if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
      echo sha256 mismatch
      exit 1
  fi

  hdiutil attach -mountpoint /Volumes/gfortran gfortran.dmg
  sudo installer -pkg /Volumes/gfortran/gfortran.pkg -target /
  otool -L /usr/local/gfortran/lib/libgfortran.3.dylib
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
  # required so that gfortran knows where to find the linking libraries.
  # export SDKROOT=/Applications/Xcode_13.3.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk
  # export SDKROOT=$(xcrun --show-sdk-path)
  type -p gfortran
fi
