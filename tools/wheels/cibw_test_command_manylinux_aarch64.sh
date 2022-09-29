set -xe

PROJECT_DIR="$1"

# python $PROJECT_DIR/tools/wheels/check_license.py
if [[ $(uname) == "Linux" || $(uname) == "Darwin" ]] ; then
    python $PROJECT_DIR/tools/openblas_support.py --check_version
fi
echo $?

# only run a reduced set of tests for cross-compiled manylinux_aarch64
python -c "import sys; import scipy.linalg; import scipy.optimize; sys.exit(not scipy.linalg.test())"
echo $?
