set -xe

python -c "import sys; import scipy; sys.exit(not scipy.test())"
