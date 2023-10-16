set -xe


NIGHTLY_FLAG=""

if [ "$#" -eq 1 ]; then
    PROJECT_DIR="$1"
elif [ "$#" -eq 2 ] && [ "$1" = "--nightly" ]; then
    NIGHTLY_FLAG="--nightly"
    PROJECT_DIR="$2"
else
    echo "Usage: $0 [--nightly] <project_dir>"
    exit 1
fi
printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt

if [ NIGHTLY_FLAG=="--nightly"];then
  python -m pip install --pre --upgrade --timeout=60 -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple scipy-openblas32
else
  python -m pip install scipy_openblas32
fi
PYTHONPATH=tools/wheels python -c "import openblas; openblas.configure_scipy_openblas()"
