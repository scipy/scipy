set -xe


if [ "$#" -eq 1 ]; then
    python -m pip insta
    python -m pip install -r requirements/openblas_requirements.txt
    PROJECT_DIR="$1"
elif [ "$#" -eq 2 ] && [ "$1" = "--nightly" ]; then
    python -m pip install --pre --upgrade --timeout=60 -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple scipy-openblas32
    PROJECT_DIR="$2"
else
    echo "Usage: $0 [--nightly] <project_dir>"
    exit 1
fi

printenv

# remove any cruft from a previous run
rm -rf build

# Update license
echo "" >> $PROJECT_DIR/LICENSE.txt
echo "----" >> $PROJECT_DIR/LICENSE.txt
echo "" >> $PROJECT_DIR/LICENSE.txt
if [[ $RUNNER_OS == "Linux" ]] ; then
    cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt
elif [[ $RUNNER_OS == "macOS" ]]; then
    cat $PROJECT_DIR/tools/wheels/LICENSE_osx.txt >> $PROJECT_DIR/LICENSE.txt
elif [[ $RUNNER_OS == "Windows" ]]; then
    cat $PROJECT_DIR/tools/wheels/LICENSE_win32.txt >> $PROJECT_DIR/LICENSE.txt
fi


# Install Openblas
echo PKG_CONFIG_PATH $PKG_CONFIG_PATH
PKG_CONFIG_PATH=$PROJECT_DIR/.openblas
rm -rf $PKG_CONFIG_PATH
mkdir -p $PKG_CONFIG_PATH
python -m pip install -r requirements/openblas_requirements.txt
python -c "import scipy_openblas64; print(scipy_openblas64.get_pkg_config())" > $PKG_CONFIG_PATH/scipy-openblas.pc
# Copy the shared objects to a path under $PKG_CONFIG_PATH, the build
# will point $LD_LIBRARY_PATH there and then auditwheel/delocate-wheel will
# pull these into the wheel. Use python to avoid windows/posix problems
python <<EOF
import os, scipy_openblas64, shutil
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), "lib")
shutil.copytree(srcdir, os.path.join("$PKG_CONFIG_PATH", "lib"))
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), ".dylibs")
if os.path.exists(srcdir):  # macosx delocate
    shutil.copytree(srcdir, os.path.join("$PKG_CONFIG_PATH", ".dylibs"))
EOF

# Install GFortran
PLATFORM=$(uname)
if [[ "$PLATFORM" == "Darwin" ]]; then
    if [[ "$(uname -m)" == "x86_64" ]]; then
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

    elif [[ $(uname -m) == "arm64" ]]; then
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
    export SDKROOT=${SDKROOT:-$(xcrun --show-sdk-path)}
elif [[ "$PLATFORM" == "Linux" ]]; then
    # nothing, use this instead of trying to detect windows with `uname`
else
    # Windows
    python -m pip install delvewheel
fi
