#!/bin/bash

# Vendors `uarray` from https://github.com/Quansight-Labs/uarray

set -o nounset
set -o errexit

REPO_URL="https://github.com/Quansight-Labs/uarray"
COMMIT_HASH="50c9cfc35e9c8f9f9b5a7499f71c34e096b3239a"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_uarray.sh`
ROOT_DIR="subprojects/uarray/uarray"

rm -rf $ROOT_DIR
mkdir -p $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
git apply ../../patches/no_setuptools_scm.patch
rm -rf src/uarray/tests
popd
mv -v $ROOT_DIR/.tmp/src/* $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/LICENSE $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp
