set -xe

WHEEL="$1"
DEST_DIR="$2"

# create a temporary directory in the destination folder and unpack the wheel
# into there
pushd $DEST_DIR
mkdir -p tmp
pushd tmp
wheel unpack $WHEEL
pushd scipy*

# To avoid DLL hell, the file name of libopenblas that's being vendored with
# the wheel has to be name-mangled. delvewheel is unable to name-mangle PYD
# containing extra data at the end of the binary, which frequently occurs when
# building with mingw.
# We therefore find each PYD in the directory structure and strip them.

for f in $(find ./scipy* -name '*.pyd'); do strip $f; done


# the libopenblas.dll is placed into this directory in the cibw_before_build
# script.
cp /c/opt/openblas/openblas_dll/*.dll ..
# Shared library for special function error handling must be copied out of
# the wheel and then reincluded with delvewheel. One cannot currently embed
# a shared library in a wheel directly on Windows.
cp ./scipy/special/libsf_error_state.dll ..

# now repack the wheel and overwrite the original
wheel pack .
mv -fv *.whl $WHEEL

cd $DEST_DIR
delvewheel repair --add-path tmp -w $DEST_DIR $WHEEL

rm -rf tmp
