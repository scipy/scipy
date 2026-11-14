#!/bin/bash

# Vendors HiGHS from https://github.com/ERGO-Code/HiGHS

set -o nounset
set -o errexit

REPO_URL="https://github.com/ERGO-Code/HiGHS"
COMMIT_HASH="f3f113d230429a8c0014b12ab0dd1a3847e58c55"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_highs.sh`
ROOT_DIR="subprojects/highs"

# start from a fresh dir
rm -rf $ROOT_DIR
mkdir $ROOT_DIR
# grab upstream into a temporary dir
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
popd
# vendor only what is needed for the Meson build, plus licenses
mv -v $ROOT_DIR/.tmp/LICENSE.txt $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/THIRD_PARTY_NOTICES.md $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/README.md $ROOT_DIR/
# METIS (Apache-2.0) is not used by the Meson build
rm -rf $ROOT_DIR/.tmp/extern/metis
mv -v $ROOT_DIR/.tmp/meson.build $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/meson_options.txt $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/highs $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/extern $ROOT_DIR/
mkdir $ROOT_DIR/highspy
mv -v $ROOT_DIR/.tmp/highspy/highs_bindings.cpp $ROOT_DIR/highspy/
# tidy up
rm -rf $ROOT_DIR/.tmp
