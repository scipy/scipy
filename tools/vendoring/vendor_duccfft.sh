# /// conda-script
# channels = ["https://prefix.dev/conda-forge"]
# entrypoint = "brush -x ${SCRIPT}"
#
# [dependencies]
# brush = "*"
# uutils-coreutils = "*"
# /// end-conda-script

# Vendors duccfft from https://gitlab.mpcdf.mpg.de/mtr/ducc.git
# Can be run via `pixi run --script tools/vendoring/vendor_duccfft.sh`
# Must be run from the repo root

set -o nounset
set -o errexit

REPO_URL="https://gitlab.mpcdf.mpg.de/mtr/ducc.git"
COMMIT_HASH="64f42ba531f609ba7029c82207a063b17f9d5275"

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
