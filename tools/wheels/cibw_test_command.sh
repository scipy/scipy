set -xe

PROJECT_DIR="$1"

# python $PROJECT_DIR/tools/wheels/check_license.py
if [[ $UNAME == "Linux" || $UNAME == "Darwin" ]] ; then
    python $PROJECT_DIR/tools/openblas_support.py --check_version
fi

python -c "import scipy; scipy.test()"
