"""
Script for checking if all the test files are installed after building.

Examples::

    $ python check_installation.py install_directory_name

        install_directory_name:
            the relative path to the directory where SciPy is installed after
            building and running `meson install`.

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


def main(install_dir):
    INSTALLED_DIR = os.path.join(ROOT_DIR, install_dir)
    if not os.path.exists(INSTALLED_DIR):
        raise ValueError(f"Provided install dir {INSTALLED_DIR} does not exist")

    scipy_test_files = get_test_files(SCIPY_DIR)
    installed_test_files = get_test_files(INSTALLED_DIR)

    # Check test files detected in repo are installed
    for test_file in scipy_test_files.keys():
        if not test_file in installed_test_files.keys():
            raise Exception("%s is not installed" % scipy_test_files[test_file])

    print("----------- All the test files were installed --------------")

    scipy_pyi_files = get_pyi_files(SCIPY_DIR)
    installed_pyi_files = get_pyi_files(INSTALLED_DIR)

    # Check *.pyi files detected in repo are installed
    for pyi_file in scipy_pyi_files.keys():
        if pyi_file not in installed_pyi_files.keys():
            raise Exception("%s is not installed" % scipy_pyi_files[pyi_file])

    print("----------- All the .pyi files were installed --------------")


def get_suffix_path(current_path, levels=1):
    current_new = current_path
    for i in range(levels + 1):
        current_new = os.path.dirname(current_new)

    return os.path.relpath(current_path, current_new)


def get_test_files(dir):
    test_files = dict()
    for path in glob.glob(f'{dir}/**/test_*.py', recursive=True):
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
    if not len(sys.argv) == 2:
        raise ValueError("Incorrect number of input arguments, need "
                         "check_installation.py relpath/to/installed/scipy")

    install_dir = sys.argv[1]
    main(install_dir)
