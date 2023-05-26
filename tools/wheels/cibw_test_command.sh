set -xe

PROJECT_DIR="$1"

# python $PROJECT_DIR/tools/wheels/check_license.py
if [[ $(uname) == "Linux" || $(uname) == "Darwin" ]] ; then
    python $PROJECT_DIR/tools/openblas_support.py --check_version
fi
echo $?

python -c "import sys; import scipy; sys.exit(not scipy.test())"
echo $?
