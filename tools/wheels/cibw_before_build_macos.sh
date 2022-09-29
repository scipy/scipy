set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")
echo $PLATFORM

# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_osx.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
basedir=$(python tools/openblas_support.py)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include

if [[ $RUNNER_OS == "macOS" && $PLATFORM == "macosx-arm64" ]]; then
  # this version of openblas has the pkg-config file included. The version
  # obtained from the openblas_support.py doesn't.
  # Problems were experienced with meson->cmake detection of openblas when
  # trying to cross compile.
  curl -L https://anaconda.org/multibuild-wheels-staging/openblas-libs/v0.3.20-140-gbfd9c1b5/download/openblas-v0.3.20-140-gbfd9c1b5-macosx_11_0_arm64-gf_f26990f.tar.gz -o openblas.tar.gz
  sudo tar -xv -C / -f openblas.tar.gz

  sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
  sudo chown -R $USER /opt/arm64-builds
  cp -r $basedir/lib/* /opt/arm64-builds/lib
  cp $basedir/include/* /opt/arm64-builds/include
fi

#########################################################################################
# Install GFortran

if [[ $PLATFORM == "macosx-x86_64" ]]; then
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

# arm64 stuff from gfortran_utils
if [[ $PLATFORM == "macosx-arm64" ]]; then
    source $PROJECT_DIR/tools/wheels/gfortran_utils.sh
    export MACOSX_DEPLOYMENT_TARGET=11.0

    # The install script requires the PLAT variable in order to set
    # the FC variable
    export PLAT=arm64
    install_arm64_cross_gfortran
    export FC=$FC_ARM64
    export PATH=$FC_LOC:$PATH
    # force a dynamic link, there may be a more elegant way of doing this.
    rm /opt/arm64-builds/lib/*.a

    # required so that gfortran knows where to find the linking libraries.
    export SDKROOT=/Applications/Xcode_13.2.1.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk
    #    export SDKROOT=$(xcrun --show-sdk-path)
    export FFLAGS=" -arch arm64 $FFLAGS"
    export LDFLAGS=" $FC_ARM64_LDFLAGS $LDFLAGS -L/opt/arm64-builds/lib -arch arm64"
    sudo ln -s $FC $FC_LOC/gfortran
    echo $(type -p gfortran)

    # having a test fortran program has helped in debugging problems with the
    # compiler environment.
    $FC $FFLAGS $PROJECT_DIR/tools/wheels/test.f $LDFLAGS
    ls -al *.out
    otool -L a.out
fi
