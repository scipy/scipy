#!/bin/bash

# Vendors duccfft from https://gitlab.mpcdf.mpg.de/mtr/ducc.git

set -o nounset
set -o errexit

REPO_URL="https://gitlab.mpcdf.mpg.de/mtr/ducc.git"
COMMIT_HASH="837c60ac4b28801cb6b17113495d4bf4cbd82968"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_pyprima.sh`
ROOT_DIR="subprojects/duccfft/ducc0"

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
popd
mv -v $ROOT_DIR/.tmp/src/ducc0/bindings/ $ROOT_DIR/bindings/
mv -v $ROOT_DIR/.tmp/src/ducc0/fft/ $ROOT_DIR/fft/
mv -v $ROOT_DIR/.tmp/src/ducc0/infra/ $ROOT_DIR/infra/
mv -v $ROOT_DIR/.tmp/src/ducc0/math/ $ROOT_DIR/math/
rm -rf $ROOT_DIR/.tmp
