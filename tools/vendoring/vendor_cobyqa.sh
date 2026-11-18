#!/bin/bash

# Vendors cobyqa from https://github.com/cobyqa/cobyqa

set -o nounset
set -o errexit

REPO_URL="https://github.com/cobyqa/cobyqa"
COMMIT_HASH="003747e754354a77bbd0ceb16c62483d2c409183" # v1.1.4

# XXX: run this from the repo top level like `./tools/vendoring/vendor_cobyqa.sh`
ROOT_DIR="subprojects/cobyqa"

# start from a fresh dir
rm -rf $ROOT_DIR
mkdir $ROOT_DIR
# grab upstream into a temporary dir
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
rm -rf cobyqa/tests
rm -rf cobyqa/utils/tests
popd
# vendor the upstream Meson build file, the package, and the license
mv -v $ROOT_DIR/.tmp/meson.build $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/LICENSE $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/cobyqa/ $ROOT_DIR/cobyqa/
# tidy up
rm -rf $ROOT_DIR/.tmp
