# /// conda-script
# channels = ["https://prefix.dev/conda-forge"]
# entrypoint = "brush -x ${SCRIPT}"
#
# [dependencies]
# brush = "*"
# uutils-coreutils = "*"
# /// end-conda-script

# Vendors biteopt from https://github.com/avaneev/biteopt
# Can be run via `pixi run --script tools/vendoring/vendor_biteopt.sh`
# Must be run from the repo root

set -o nounset
set -o errexit

REPO_URL="https://github.com/avaneev/biteopt"
COMMIT_HASH="9ccb2352443d8472a4675b6e6f92bff2adaeaea7"

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
