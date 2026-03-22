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
git grep -l "SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later" | xargs tar cf ducc_bsd.tar
popd
tar xf $ROOT_DIR/.tmp/ducc_bsd.tar
rm -rf $ROOT_DIR/.tmp
