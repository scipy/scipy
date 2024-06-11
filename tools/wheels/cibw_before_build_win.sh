set -xe

PROJECT_DIR="$1"

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_win32.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
python -m pip install -r requirements/openblas.txt
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PROJECT_DIR/scipy-openblas.pc

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel
