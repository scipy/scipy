#!/usr/bin/env python
import os
import sys
import subprocess
from argparse import ArgumentParser


CONFIG = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    'lint.ini',
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
        ['git', 'diff', '--name-only', '-z', sha, '--',
         '*.py', '*.pyx', '*.pxd', '*.pxi'],
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    res.check_returncode()
    return [f for f in res.stdout.split('\0') if f]


def run_flake8(files):
    if not files:
        return 0, ""
    res = subprocess.run(
        ['flake8', '--config', CONFIG] + files,
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    return res.returncode, res.stdout


def main():
    parser = ArgumentParser()
    # In Python 3.9, can use: argparse.BooleanOptionalAction
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

    rc, errors = run_flake8(files)

    if errors:
        print(errors)

    sys.exit(rc)


if __name__ == '__main__':
    main()
