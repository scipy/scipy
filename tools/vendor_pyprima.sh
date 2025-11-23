#!/bin/bash

# Vendors pyprima from https://github.com/libprima/prima

set -o nounset
set -o errexit

REPO_URL="https://github.com/libprima/prima"
COMMIT_HASH="6fbf15d9b4aac20d2e3ee5b3fcb41f0ca17c4330"

# XXX: run this from the repo top level like `./tools/vendor_pyprima.sh`
ROOT_DIR="scipy/_lib/pyprima"

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
rm -rf pyprima/profiles_vs_matlab
rm -rf pyprima/profiles_vs_python_bindings
popd
mv -v $ROOT_DIR/.tmp/LICENCE.txt $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/pyprima/ $ROOT_DIR/pyprima/
rm -rf $ROOT_DIR/.tmp
