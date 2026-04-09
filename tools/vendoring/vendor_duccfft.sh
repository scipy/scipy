#!/bin/bash

# Vendors duccfft from https://gitlab.mpcdf.mpg.de/mtr/ducc.git

set -o nounset
set -o errexit

REPO_URL="https://gitlab.mpcdf.mpg.de/mtr/ducc.git"
COMMIT_HASH="8e305676771659bbb562c4c7d0c59b921ba46da8"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_duccfft.sh`
ROOT_DIR="subprojects/duccfft/ducc0"

# start from a fresh dir
rm -rf $ROOT_DIR
# create needed directories
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
mkdir $ROOT_DIR/.tmpBSD
# grab upstream into a temporary dir
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
# extract license-compatible code
git grep -l "SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later" | xargs tar cf ducc_bsd.tar
popd
tar xf $ROOT_DIR/.tmp/ducc_bsd.tar -C $ROOT_DIR/.tmpBSD
# vendor code into the final location
mv -v $ROOT_DIR/.tmpBSD/src/ducc0/* $ROOT_DIR/
# tidy up
rm -rf $ROOT_DIR/.tmp
rm -rf $ROOT_DIR/.tmpBSD
