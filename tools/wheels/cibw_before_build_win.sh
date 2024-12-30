set -xe

PROJECT_DIR="$1"

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_win32.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
python -m pip install -r requirements/openblas.txt
python -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > $PROJECT_DIR/scipy-openblas.pc

# TODO: delete along with enabling build isolation by unsetting
# CIBW_BUILD_FRONTEND when scipy is buildable under free-threaded
# python with a released version of cython
FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    python -m pip install -U --pre pip
    python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy cython
    python -m pip install git+https://github.com/serge-sans-paille/pythran
    python -m pip install ninja meson-python pybind11
fi

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel wheel
