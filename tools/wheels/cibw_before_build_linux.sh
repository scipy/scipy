set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
basedir=$(python tools/openblas_support.py)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
if [[ $RUNNER_OS == "macOS" && $PLATFORM == "macosx-arm64" ]]; then
    sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
    sudo chown -R $USER /opt/arm64-builds
    cp -r $basedir/lib/* /opt/arm64-builds/lib
    cp $basedir/include/* /opt/arm64-builds/include
fi
