set -xe

PROJECT_DIR="$1"

python -c "import sys; import scipy; sys.exit(not scipy.test('full', extra_argv=['-vvv']))"

python $PROJECT_DIR/tools/wheels/check_license.py
if [[ $UNAME == "Linux" || $UNAME == "Darwin" ]] ; then
    python $PROJECT_DIR/tools/openblas_support.py --check_version
fi
