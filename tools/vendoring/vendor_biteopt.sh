#!/bin/bash

# Vendors biteopt from https://github.com/avaneev/biteopt

set -o nounset
set -o errexit

REPO_URL="https://github.com/avaneev/biteopt"
COMMIT_HASH="9ccb2352443d8472a4675b6e6f92bff2adaeaea7"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_biteopt.sh`
ROOT_DIR="subprojects/biteopt/biteopt"

# start from a fresh dir
rm -rf $ROOT_DIR
mkdir $ROOT_DIR
# grab upstream into a temporary dir
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
# vendor the include closure of biteopt.h (header-only dependency) and the license
mv -v biteopt.h spheropt.h mbopt.h nmsopt.h biteaux.h LICENSE ..
popd
# tidy up
rm -rf $ROOT_DIR/.tmp
