#!/usr/bin/env python
import os
import sys
import subprocess
import packaging.version
from argparse import ArgumentParser


CONFIG = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    'lint.toml',
)


def rev_list(branch, num_commits):
    """List commits in reverse chronological order.

    Only the first `num_commits` are shown.

    """
    res = subprocess.run(
        [
            'git',
            'rev-list',
            '--max-count',
            f'{num_commits}',
            '--first-parent',
            branch
        ],
        stdout=subprocess.PIPE,
        encoding='utf-8',
    )
    res.check_returncode()
    return res.stdout.rstrip('\n').split('\n')


def find_branch_point(branch):
    """Find when the current branch split off from the given branch.

    It is based off of this Stackoverflow post:

    https://stackoverflow.com/questions/1527234/finding-a-branch-point-with-git#4991675

    """
    branch_commits = rev_list('HEAD', 1000)
    main_commits = set(rev_list(branch, 1000))
    for branch_commit in branch_commits:
        if branch_commit in main_commits:
            return branch_commit

    # If a branch split off over 1000 commits ago we will fail to find
    # the ancestor.
    raise RuntimeError(
        'Failed to find a common ancestor in the last 1000 commits'
    )


def diff_files(sha):
    """Find the diff since the given SHA."""
    res = subprocess.run(
        ['git', 'diff', '--name-only', '--diff-filter=ACMR', '-z', sha, '--',
         '*.py', '*.pyx', '*.pxd', '*.pxi'],
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    res.check_returncode()
    return [f for f in res.stdout.split('\0') if f]


def run_ruff(files, fix):
    if not files:
        return 0, ""
    args = ['--fix', '--exit-non-zero-on-fix'] if fix else []
    res = subprocess.run(
        ['ruff', 'check', f'--config={CONFIG}'] + args + list(files),
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    return res.returncode, res.stdout


def run_cython_lint(files):
    if not files:
        return 0, ""
    res = subprocess.run(
        ['cython-lint', '--no-pycodestyle'] + list(files),
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    return res.returncode, res.stdout


def check_ruff_version():
    min_version = packaging.version.parse('0.0.292')
    res = subprocess.run(
        ['ruff', '--version'],
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    version = res.stdout.replace('ruff ', '')
    if packaging.version.parse(version) < min_version:
        raise RuntimeError("Linting requires `ruff>=0.0.292`. Please upgrade `ruff`.")


def main():
    check_ruff_version()
    parser = ArgumentParser(description="Also see `pre-commit-hook.py` which "
                                        "lints all files staged in git.")
    # In Python 3.9, can use: argparse.BooleanOptionalAction
    parser.add_argument("--fix", action='store_true',
                        help='Attempt to fix linting violations')
    parser.add_argument("--diff-against", dest='branch',
                        type=str, default=None,
                        help="Diff against "
                             "this branch and lint modified files. Use either "
                             "`--diff-against` or `--files`, but not both.")
    parser.add_argument("--files", nargs='*',
                        help="Lint these files or directories; "
                             "use **/*.py to lint all files")

    args = parser.parse_args()

    if not ((args.files is None) ^ (args.branch is None)):
        print('Specify either `--diff-against` or `--files`. Aborting.')
        sys.exit(1)

    if args.branch:
        branch_point = find_branch_point(args.branch)
        files = diff_files(branch_point)
    else:
        files = args.files

    cython_exts = ('.pyx', '.pxd', '.pxi')
    cython_files = {f for f in files if any(f.endswith(ext) for ext in cython_exts)}
    other_files = set(files) - cython_files

    rc_cy, errors = run_cython_lint(cython_files)
    if errors:
        print(errors)

    rc, errors = run_ruff(other_files, fix=args.fix)
    if errors:
        print(errors)

    if rc == 0 and rc_cy != 0:
        rc = rc_cy

    sys.exit(rc)


if __name__ == '__main__':
    main()
