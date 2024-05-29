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

# Install Openblas
python -m pip install -r requirements/openblas.txt
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PROJECT_DIR/scipy-openblas.pc
