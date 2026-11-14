#!/bin/bash

# Vendors SuperLU from https://github.com/xiaoyeli/superlu

set -o nounset
set -o errexit

REPO_URL="https://github.com/xiaoyeli/superlu"
COMMIT_HASH="a8d2543611ba38b129b619751ca2a7284e624a07"

# XXX: run this from the repo top level like `./tools/vendoring/vendor_superlu.sh`
ROOT_DIR="subprojects/superlu/superlu"

# start from a fresh dir
rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/SRC
# grab upstream into a temporary dir
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
# `mc64ad.c` is not license-compatible with SciPy, so don't vendor it
rm SRC/mc64ad.c
popd
# vendor the license, dropping the "External software" section: it only
# covers the MC64 routine in `mc64ad.c`, which we don't vendor
python -I - $ROOT_DIR/.tmp/License.txt > $ROOT_DIR/License.txt <<'EOF'
import sys
text = open(sys.argv[1]).read()
print(text.split("-----------------\nExternal software")[0].rstrip("\n"))
EOF
# vendor the README, library sources, and config header template
mv -v $ROOT_DIR/.tmp/README $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/SRC/*.c $ROOT_DIR/.tmp/SRC/*.h $ROOT_DIR/SRC/
mv -v $ROOT_DIR/.tmp/SRC/superlu_config.h.in $ROOT_DIR/SRC/
# tidy up
rm -rf $ROOT_DIR/.tmp
