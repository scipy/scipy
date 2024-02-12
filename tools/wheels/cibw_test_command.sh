set -xe

PROJECT_DIR="$1"

# python $PROJECT_DIR/tools/wheels/check_license.py

python -c "import sys; import scipy; sys.exit(not scipy.test())"
echo $?
