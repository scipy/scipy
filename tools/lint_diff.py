#!/usr/bin/env python
import os
import sys
import subprocess
from argparse import ArgumentParser

CONFIG = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
    'lint_diff.ini',
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
        'Failed to find a common ancestor in the last 1000 commits')


def find_diff(sha, files=None):
    """Find the diff since the given sha."""
    if files:
        for file_or_dir in files:
            msg = f"{file_or_dir} doesn't exist. Please provide a valid path."
            assert os.path.exists(file_or_dir), msg
    else:
        files = ['*.py']
    res = subprocess.run(
        ['git', 'diff', '--unified=0', sha, '--'] + files,
        stdout=subprocess.PIPE,
        encoding='utf-8'
    )
    res.check_returncode()
    return res.stdout


def run_flake8(diff):
    """Run flake8 on the given diff."""
    res = subprocess.run(
        ['flake8', '--diff', '--config', CONFIG],
        input=diff,
        stdout=subprocess.PIPE,
        encoding='utf-8',
    )
    return res.returncode, res.stdout


def main():
    parser = ArgumentParser()
    parser.add_argument("--branch", type=str, default='main',
                        help="The branch to diff against")
    parser.add_argument("--files", type=str, nargs='+', default=None,
                        help="The files or directories to diff against")
    args = parser.parse_args()

    branch_point = find_branch_point(args.branch)
    diff = find_diff(branch_point, args.files)
    rc, errors = run_flake8(diff)
    if errors:
        print(errors)
    else:
        print("No lint errors found.")
    sys.exit(rc)


if __name__ == '__main__':
    main()
