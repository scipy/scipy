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

  source tools/wheels/gfortran_utils.sh
  # for this compiler to work SDKROOT needs to be set, along the lines of
  # export SDKROOT=/Applications/Xcode_11.7.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk
  install_gfortran

  type -p gfortran
  gfortran --version
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
  # export SDKROOT=/Applications/Xcode_13.2.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk
  # export SDKROOT=$(xcrun --show-sdk-path)
  type -p gfortran
fi
