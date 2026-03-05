#!/bin/bash

# Vendors `uarray` from https://github.com/Quansight-Labs/uarray

set -o nounset
set -o errexit

REPO_URL="https://github.com/Quansight-Labs/uarray"
COMMIT_HASH="07f399bd7a3566bd5e3cba7f600463bbf4672e03"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_uarray.sh`
ROOT_DIR="subprojects/uarray"

rm -rf $ROOT_DIR
mkdir -p $ROOT_DIR/.tmp1
mkdir -p $ROOT_DIR/.tmp2

# clone source, build and extract sdist
git clone $REPO_URL $ROOT_DIR/.tmp1
pushd $ROOT_DIR/.tmp1
git checkout $COMMIT_HASH
uv build --sdist
SDIST="$(ls dist/*.tar.gz)"
tar -xzf "$SDIST" -C ../.tmp2 --strip-components=1
popd
rm -rf $ROOT_DIR/.tmp1

# grab the files we want from the sdist
mkdir -p $ROOT_DIR/src
mv -v $ROOT_DIR/.tmp2/src/* $ROOT_DIR/src/
mv -v $ROOT_DIR/.tmp2/LICENSE $ROOT_DIR/
mv -v $ROOT_DIR/.tmp2/meson.build $ROOT_DIR/
mv -v $ROOT_DIR/.tmp2/pyproject.toml $ROOT_DIR/
mkdir -p $ROOT_DIR/tools
mv -v $ROOT_DIR/.tmp2/tools/gitversion.py $ROOT_DIR/tools
rm -rf $ROOT_DIR/.tmp2

cd $ROOT_DIR
echo 'This directory is populated by `tools/vendoring/vendor_uarray.sh`.' > README.md
