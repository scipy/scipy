#!/bin/bash

# Vendors `packaging.version` from https://github.com/pypa/packaging

set -o nounset
set -o errexit

REPO_URL="https://github.com/pypa/packaging"
COMMIT_HASH="f58537628042c7f29780b9d33f31597e7fc9d664"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_packaging_version.sh`
ROOT_DIR="scipy/_external/packaging_version/src"

rm -rf $ROOT_DIR
mkdir -p $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
popd
mv -v $ROOT_DIR/.tmp/src/packaging/version.py $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/src/packaging/_structures.py $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/LICENSE.BSD $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp
