#!/bin/bash

# Vendors qhull from https://github.com/qhull/qhull

set -o nounset
set -o errexit

REPO_URL="https://github.com/qhull/qhull"
# commit for v8.1-alpha6 of qhull
COMMIT_HASH="d1c2fc0caa5f644f3a0f220290d4a868c68ed4f6"

# XXX: run this from the repo top level like `./tools/vendor_qhull.sh`
ROOT_DIR="subprojects/qhull_r/libqhull_r"

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
pushd src/libqhull_r/
rm *.htm
rm *.pro
rm *.def
rm Makefile
popd # $ROOT_DIR/.tmp
popd
mv -v $ROOT_DIR/.tmp/COPYING.txt $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/Announce.txt $ROOT_DIR/
cp -v $ROOT_DIR/.tmp/src/libqhull_r/* $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp
