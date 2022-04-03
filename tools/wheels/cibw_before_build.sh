set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

# Install Openblas
basedir=$(python tools/openblas_support.py)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
