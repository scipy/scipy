set -ex

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: delete when numpy is buildable under free-threaded python
    # with a released version of cython
    python -m pip install -U --pre pip
    python -m pip install git+https://github.com/cython/cython
    python -m pip install -i https://pypi.anaconda.org/scientific-python-nightly-wheels/simple numpy
fi
