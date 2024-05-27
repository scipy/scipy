"""
Script for checking if all the test files are installed after building.

Examples::

    $ python check_installation.py install_directory_name

        install_directory_name:
            the relative path from the root of the repo to the directory where
            SciPy is installed (for dev.py usually "build-install")

Notes
=====

The script will stop on encountering the first missing file in the install dir,
it will not give a full listing. This should be okay, because the script is
meant for use in CI so it's not like many files will be missing at once.

"""

import os
import glob
import sys


CUR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))
ROOT_DIR = os.path.dirname(CUR_DIR)
SCIPY_DIR = os.path.join(ROOT_DIR, 'scipy')


# Files whose installation path will be different from original one
changed_installed_path = {
    'scipy/_build_utils/tests/test_scipy_version.py':
        'scipy/_lib/tests/test_scipy_version.py'
}

# We do not want the following tests to be checked.
# Note: only 3 subdirs are listed here, due to how `get_test_files` is
# implemented.
# If this list gets too annoying, we should implement excluding directories
# rather than (or in addition to) files.
exception_list_test_files = [
    "_lib/array_api_compat/tests/test_all.py",
    "_lib/array_api_compat/tests/test_array_namespace.py",
    "_lib/array_api_compat/tests/test_common.py",
    "_lib/array_api_compat/tests/test_isdtype.py",
    "_lib/array_api_compat/tests/test_vendoring.py",
    "_lib/array_api_compat/tests/test_array_namespace.py",
    "cobyqa/cobyqa/tests/test_main.py",
    "cobyqa/cobyqa/tests/test_models.py",
    "cobyqa/cobyqa/tests/test_problem.py",
    "cobyqa/utils/tests/test_exceptions.py",
    "cobyqa/utils/tests/test_math.py",
]


def main(install_dir, no_tests):
    INSTALLED_DIR = os.path.join(ROOT_DIR, install_dir)
    if not os.path.exists(INSTALLED_DIR):
        raise ValueError(f"Provided install dir {INSTALLED_DIR} does not exist")

    scipy_test_files = get_test_files(SCIPY_DIR)
    scipy_test_extension_modules = get_test_files(INSTALLED_DIR, "so")
    installed_test_files = get_test_files(INSTALLED_DIR)

    if no_tests:
        if len(scipy_test_extension_modules) > 0:
            raise Exception(f"{scipy_test_extension_modules.values()} "
                            "should not be installed but "
                            "are found in the installation directory.")
    else:
        if len(scipy_test_extension_modules) == 0:
            raise Exception("Test for extension modules should be "
                            "installed but are not found in the "
                            "installation directory.")

    # Check test files detected in repo are installed
    for test_file in scipy_test_files.keys():
        if test_file in exception_list_test_files:
            continue

        if no_tests:
            if test_file in installed_test_files:
                raise Exception(f"{test_file} should not be installed but "
                        "is found in the installation directory.")
            continue

        if test_file not in installed_test_files.keys():
            raise Exception(f"{scipy_test_files[test_file]} is not installed; "
                            f"either install it or add `{test_file}` to the "
                            "exception list in `tools/check_installation.py`")

    if no_tests:
        print("----------- No test files were installed --------------")
    else:
        print("----------- All the test files were installed --------------")

    scipy_pyi_files = get_pyi_files(SCIPY_DIR)
    installed_pyi_files = get_pyi_files(INSTALLED_DIR)

    # Check *.pyi files detected in repo are installed
    for pyi_file in scipy_pyi_files.keys():
        if pyi_file not in installed_pyi_files.keys():
            if no_tests and "test" in scipy_pyi_files[pyi_file]:
                continue
            raise Exception("%s is not installed" % scipy_pyi_files[pyi_file])

    print("----------- All the necessary .pyi files were installed --------------")


def get_suffix_path(current_path, levels=1):
    current_new = current_path
    for i in range(levels + 1):
        current_new = os.path.dirname(current_new)

    return os.path.relpath(current_path, current_new)


def get_test_files(dir, ext="py"):
    test_files = dict()
    underscore = "_" if ext == "so" else ""
    for path in glob.glob(f'{dir}/**/{underscore}test_*.{ext}', recursive=True):
        suffix_path = get_suffix_path(path, 3)
        suffix_path = changed_installed_path.get(suffix_path, suffix_path)
        if "highspy" not in suffix_path:
            test_files[suffix_path] = path

    return test_files


def get_pyi_files(dir):
    pyi_files = dict()
    for path in glob.glob(f'{dir}/**/*.pyi', recursive=True):
        suffix_path = get_suffix_path(path, 2)
        pyi_files[suffix_path] = path

    return pyi_files


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise ValueError("Incorrect number of input arguments, need "
                         "check_installation.py relpath/to/installed/scipy")

    install_dir = sys.argv[1]
    no_tests = False
    if len(sys.argv) == 3:
        no_tests = sys.argv[2] == "--no-tests"
    main(install_dir, no_tests)
