#!/bin/bash

# Vendors cobyqa from https://github.com/cobyqa/cobyqa

set -o nounset
set -o errexit

REPO_URL="https://github.com/cobyqa/cobyqa"
COMMIT_HASH="af497d2a89dc9db233053615da20a1c687ec4c1f"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_cobyqa.sh`
ROOT_DIR="subprojects/cobyqa/cobyqa"

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
rm -rf cobyqa/tests
rm -rf cobyqa/utils/tests
popd
mv -v $ROOT_DIR/.tmp/LICENSE $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/cobyqa/ $ROOT_DIR/cobyqa/
rm -rf $ROOT_DIR/.tmp
