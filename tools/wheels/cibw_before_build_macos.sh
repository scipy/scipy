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
fi
