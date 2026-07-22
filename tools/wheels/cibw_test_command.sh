set -xe

FREE_THREADED_BUILD="$(python -c"import sysconfig; print(bool(sysconfig.get_config_var('Py_GIL_DISABLED')))")"
if [[ $FREE_THREADED_BUILD == "True" ]]; then
    # Manually check that importing SciPy does not re-enable the GIL.
    # In principle the tests should catch this but it seems harmless to leave it
    # here as a final sanity check before uploading broken wheels
    if [[ $(python -c "import scipy.stats" 2>&1) == *"The global interpreter lock (GIL) has been enabled"* ]]; then
        echo "Error: Importing SciPy re-enables the GIL in the free-threaded build"
        exit 1
    fi

fi

python -c "import sys; import scipy; sys.exit(not scipy.test())"
