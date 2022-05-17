import os
import sys
import subprocess
from pathlib import Path
import argparse
import textwrap


repo_root_dir = Path(__file__).resolve().parent.parent


def check_submodules():
    """
    Verify that the submodules are checked out and clean.

    Use `git submodule update --init` if this check fails.
    """
    if not os.path.exists(repo_root_dir / '.git'):
        return
    with open(repo_root_dir / '.gitmodules') as f:
        for line in f:
            if 'path' in line:
                p = line.split('=')[-1].strip()
                if not os.path.exists(p):
                    raise ValueError('Submodule %s missing' % p)

    proc = subprocess.Popen(['git', 'submodule', 'status'],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    for line in status.splitlines():
        if line.startswith('-') or line.startswith('+'):
            raise ValueError('Submodule not clean: %s' % line)


def check_not_dirty():
    proc = subprocess.Popen(['git', 'status'], stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    if "nothing to commit, working tree clean" not in status:
        print(f"{status}\n")
        raise RuntimeError("Git repo dirty, so cannot generate sdist for PyPI")


def concat_license_files():
    """
    Merge LICENSE.txt and LICENSES_bundled.txt for sdist creation

    Done this way to keep LICENSE.txt in repo as exact BSD 3-clause (see
    NumPy gh-13447).  This makes GitHub state correctly how SciPy is licensed.
    """
    file1 = repo_root_dir / 'LICENSE.txt'
    file2 = repo_root_dir / 'LICENSES_bundled.txt'

    # Concatenate files
    with open(file1, 'a') as f1:
        with open(file2, 'r') as f2:
            bundled_text = f2.read()
            f1.write('\n\n')
            f1.write(bundled_text)

    # remove LICENSES_bundled.txt, so we don't add it to the sdist
    os.remove(file2)


numpy_reqs_toml = """
    # NumPy dependencies - to update these, sync from
    # https://github.com/scipy/oldest-supported-numpy/, and then
    # update minimum version to match our install_requires min version
    # ----------------------------------------------------------------

    # numpy 1.19 was the first minor release to provide aarch64 wheels, but
    # wheels require fixes contained in numpy 1.19.2
    "numpy==1.19.2; python_version=='3.8' and platform_machine=='aarch64' and platform_python_implementation != 'PyPy'",
    # aarch64 for py39 is covered by default requirement below

    # arm64 on Darwin supports Python 3.8 and above requires numpy>=1.20.0
    "numpy==1.20.0; python_version=='3.8' and platform_machine=='arm64' and platform_system=='Darwin'",
    "numpy==1.20.0; python_version=='3.9' and platform_machine=='arm64' and platform_system=='Darwin'",

    # default numpy requirements
    "numpy==1.18.5; python_version=='3.8' and (platform_machine!='arm64' or platform_system!='Darwin') and platform_machine!='aarch64' and platform_python_implementation != 'PyPy'",
    "numpy==1.19.3; python_version=='3.9' and (platform_machine!='arm64' or platform_system!='Darwin') and platform_python_implementation != 'PyPy'",
    "numpy==1.21.4; python_version=='3.10' and platform_python_implementation != 'PyPy'",

    # For Python versions which aren't yet officially supported,
    # we specify an unpinned NumPy which allows source distributions
    # to be used and allows wheels to be used as soon as they
    # become available.
    "numpy; python_version>='3.11'",
    "numpy; python_version>='3.8' and platform_python_implementation=='PyPy'",
"""


def modify_pyproject_toml():
    pyproject_toml = repo_root_dir / "pyproject.toml"
    with open(pyproject_toml, 'r') as f1:
        content = f1.readlines()

    # Find unpinned numpy requirements (build and runtime)
    idx_numpyreq = []
    for idx, line in enumerate(content):
        if "numpy>=" in line:
            idx_numpyreq.append(idx)

    if len(idx_numpyreq) != 2:
        raise RuntimeError("Expected 2 lines with `numpy>=` in pyproject.toml"
                           f", found: {len(idx_numpyreq)}")

    # If we're building an sdist, we can't change the runtime req. - so only
    # update the build req. here:
    content[idx_numpyreq[0]] = numpy_reqs_toml
    with open(pyproject_toml, 'w') as f1:
        f1.write("".join(content))


def get_current_commit_id():
    git_dir = os.path.join(repo_root_dir, ".git")
    proc = subprocess.Popen(['git', '--git-dir', git_dir,
                             'rev-parse', '--verify', 'HEAD'],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    if status.endswith('\n'):
        commit_id = status[:-1]

    assert len(commit_id) == 40
    return commit_id


def commit_changed_files():
    proc = subprocess.Popen(['git', 'commit', '-a', '-m',
                             '"Commit for PyPI sdist (by build-pypi-artifacts.py)"'],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    if not ("2 files changed" in status and "rename LICENSES" in status):
        print(f"{status}\n")
        raise RuntimeError("Expected changes to 3 files: pyproject.toml, "
                           "LICENSE.txt, and LICENSES_bundled.txt")


def reset_to_previous_commit(commit_id):
    proc = subprocess.Popen(['git', 'reset', '--hard', commit_id],
                            stdout=subprocess.PIPE)
    status, _ = proc.communicate()
    status = status.decode("ascii", "replace")
    if "HEAD is now at" not in status:
        print(f"{status}\n")
        raise RuntimeError("Resetting to original commit may have gone wrong")


def generate_sdist():
    subprocess.run([sys.executable, '-m', 'build', '--sdist',
                    '--no-isolation', '--skip-dependency-check'])
    distdir = repo_root_dir / 'dist'
    if not len(list(distdir.glob('scipy-*.tar.gz'))) == 1:
        raise RuntimeError("sdist creating not succesful!")


def create_sdist_for_pypi():
    check_submodules()
    check_not_dirty()
    commit_id = get_current_commit_id()

    concat_license_files()
    modify_pyproject_toml()

    commit_changed_files()
    generate_sdist()
    reset_to_previous_commit(commit_id)


def create_wheel_for_pypi(skip_check=False):
    """
    Create a wheel with the appropriate wheel-specific license info (for
    vendored libraries) and `_distributor_init.py` file.
    """
    # TODO: adapt for what is going on in https://github.com/MacPython/scipy-wheels
    #       - grab the correct build of OpenBLAS (or keep separate?)
    #       - patch_code.sh changes need to live here (double check and
    #         document what mesonpy does with vendoring)
    #       - patch_numpy.sh isn't needed
    #       - set MACOSX_DEPLOYMENT_TARGET in environment here
    #       - unclear whether we have a use for setting PYTHONFAULTHANDLER=1
    #       - verify_init.py is likely good to incorporate here
    raise NotImplementedError('TODO')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=textwrap.indent(
            textwrap.dedent(
                """
                Tool to build SciPy release artifacts (sdist, wheels) for PyPI.

                This tool is needed for releases because we need to take
                PyPI-specific actions like amend license files and
                the numpy dependency in pyproject.toml before invoking a
                standard tool like `pypa/build`.

                Usage:

                  $ python build-pypi-artifacts.py            # build sdist, then a wheel
                  $ python build-pypi-artifacts.py --sdist    # build sdist only
                  $ python build-pypi-artifacts.py --wheel path/to/sdist

                Note that the sdist can be created in any environment, there
                is no platform or Python-version specific content that ends up
                in the sdist. The wheels on the other hand must be created in
                an environment that has the correct dependencies already
                installed. If there's a mismatch in, for example, the version
                of numpy that is installed in the (CI) environment compared to
                what's specified in `pyproject.toml`, an error will be raised.

                For testing purposes, the dependency check for wheels can be
                bypassed with `--skip-dependency-check`.
                """
            ).strip(),
            '    ',
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '--sdist',
        '-s',
        action='store_true',
        help='Build a source distribution',
    )
    parser.add_argument(
        '--wheel',
        '-w',
        type=str,
        help='Build a wheel from an sdist (provide path to sdist as argument)',
    )
    parser.add_argument(
        '--skip-dependency-check',
        '-x',
        action='store_true',
        help='Do not check that build dependencies are installed',
    )
    args = parser.parse_args()

    if args.wheel:
        # TODO: wheel build
        raise NotImplementedError('TODO')
        create_wheel_for_pypi(skip_check=args.skip_dependency_check)
    elif args.sdist:
        create_sdist_for_pypi()
    else:
        create_sdist_for_pypi()
        raise NotImplementedError('TODO')
        create_wheel_for_pypi(skip_check=args.skip_dependency_check)
