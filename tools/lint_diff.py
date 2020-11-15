#!/usr/bin/env python
import os
import sys
import subprocess
<<<<<<< HEAD
=======
from argparse import ArgumentParser
>>>>>>> 2a9e4923aa2be5cd54ccf2196fc0da32fe459e76

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


<<<<<<< HEAD
def find_branch_point():
    """Find when the current branch split off master.
=======
def find_branch_point(branch):
    """Find when the current branch split off from the given branch.
>>>>>>> 2a9e4923aa2be5cd54ccf2196fc0da32fe459e76

    It is based off of this Stackoverflow post:

    https://stackoverflow.com/questions/1527234/finding-a-branch-point-with-git#4991675

    """
    branch_commits = rev_list('HEAD', 1000)
<<<<<<< HEAD
    master_commits = set(rev_list('master', 1000))
=======
    master_commits = set(rev_list(branch, 1000))
>>>>>>> 2a9e4923aa2be5cd54ccf2196fc0da32fe459e76
    for branch_commit in branch_commits:
        if branch_commit in master_commits:
            return branch_commit

    # If a branch split off over 1000 commits ago we will fail to find
    # the ancestor.
    raise RuntimeError(
        'Failed to find a common ancestor in the last 1000 commits')


def find_diff(sha):
    """Find the diff since the given sha."""
    res = subprocess.run(
        ['git', 'diff', '--unified=0', sha, '--', '*.py'],
        stdout=subprocess.PIPE,
        encoding='utf-8',
    )
    res.check_returncode()
    return res.stdout


def run_pycodestyle(diff):
    """Run pycodestyle on the given diff."""
    res = subprocess.run(
        ['pycodestyle', '--diff', '--config', CONFIG],
        input=diff,
        stdout=subprocess.PIPE,
        encoding='utf-8',
    )
    return res.returncode, res.stdout


def main():
<<<<<<< HEAD
    branch_point = find_branch_point()
=======
    parser = ArgumentParser()
    parser.add_argument("--branch", type=str, default='master',
                        help="The branch to diff against")
    args = parser.parse_args()

    branch_point = find_branch_point(args.branch)
>>>>>>> 2a9e4923aa2be5cd54ccf2196fc0da32fe459e76
    diff = find_diff(branch_point)
    rc, errors = run_pycodestyle(diff)
    if errors:
        print(errors)
    sys.exit(rc)


if __name__ == '__main__':
    main()
