set -xe

PROJECT_DIR="$1"

# Install OpenBLAS
python -m pip install -r requirements/openblas.txt
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PROJECT_DIR/scipy-openblas.pc

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel wheel
