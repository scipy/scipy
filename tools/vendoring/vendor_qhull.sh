# /// conda-script
# channels = ["https://prefix.dev/conda-forge"]
# entrypoint = "brush -x ${SCRIPT}"
#
# [dependencies]
# brush = "*"
# uutils-coreutils = "*"
# /// end-conda-script

# Vendors qhull from https://github.com/qhull/qhull
# Can be run via `pixi run --script tools/vendoring/vendor_qhull.sh`
# Must be run from the repo root

set -o nounset
set -o errexit

REPO_URL="https://github.com/qhull/qhull"
COMMIT_HASH="613debeaea72ee66626dace9ba1a2eff11b5d37d"

ROOT_DIR="subprojects/qhull_r/libqhull_r"

rm -rf $ROOT_DIR
mkdir $ROOT_DIR
mkdir $ROOT_DIR/.tmp
git clone $REPO_URL $ROOT_DIR/.tmp
pushd $ROOT_DIR/.tmp
git checkout $COMMIT_HASH
# TODO: delete this git apply when qhull v8.2.0 is released
git apply ../../patches/poly2rc.patch
pushd src/libqhull_r/
rm *.htm
rm *.pro
rm *.def
rm Makefile
popd # $ROOT_DIR/.tmp
popd
mv -v $ROOT_DIR/.tmp/COPYING.txt $ROOT_DIR/
mv -v $ROOT_DIR/.tmp/Announce.txt $ROOT_DIR/
cp -v $ROOT_DIR/.tmp/src/libqhull_r/* $ROOT_DIR/
rm -rf $ROOT_DIR/.tmp
