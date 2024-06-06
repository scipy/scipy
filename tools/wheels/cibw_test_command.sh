set -xe

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # TODO: delete when importing numpy no longer enables the GIL
    # setting to zero ensures the GIL is disabled while running the
    # tests under free-threaded python
    export PYTHON_GIL=0
fi

python -c "import sys; import scipy; sys.exit(not scipy.test())"
