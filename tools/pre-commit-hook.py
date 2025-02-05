#!/usr/bin/env python
#
# Pre-commit linting hook.
#
# Install from root of repository with:
#
#   cp tools/pre-commit-hook.py .git/hooks/pre-commit

import subprocess
import sys
import os


# Run lint.py from the scipy source tree
linters = [
    '../../tools/lint.py',
    'tools/lint.py',
    'lint.py'  # in case pre-commit hook is run from tools dir
]

linter = [f for f in linters if os.path.exists(f)][0]

unicode_checks = [
    '../../tools/check_unicode.py',
    'tools/check_unicode.py',
    'check_unicode.py'  # in case pre-commit hook is run from tools dir
]

unicode_check = [f for f in unicode_checks if os.path.exists(f)][0]


# names of files that were staged
# add  '*.pxd', '*.pxi' once cython-lint supports it
p = subprocess.run(['git', 'diff',
                    '--cached', '--name-only', '-z',
                    '--diff-filter=ACMR',
                    '--', '*.py', '*.pyx'],
                   capture_output=True, check=True)

files = p.stdout.decode(sys.getfilesystemencoding()).split('\0')
files = [f for f in files if f]

# create a temporary copy of what would get committed, without unstaged
# modifications (e.g., only certain changes in a file may have been committed)
git_dir = os.environ.get('GIT_DIR', '.git')
work_dir = os.path.join(git_dir, '.pre-commit-work_dir')

p = subprocess.run(['git', 'write-tree'], capture_output=True, check=True)
tree_hash = p.stdout.decode('ascii').split('\n')[0]

p = subprocess.run(['git', 'commit-tree', '-p', 'HEAD',
                    tree_hash, '-m', '...'], capture_output=True, check=True)
fake_commit = p.stdout.decode('ascii').split('\n')[0]

if not os.path.isdir(work_dir):
    subprocess.run(['git', 'clone', '-qns', git_dir, work_dir])

subprocess.run(['git', 'reset', '--quiet', '--hard', 'HEAD'],
               env={}, cwd=work_dir, check=True)
subprocess.run(['git', 'checkout', '-q', fake_commit],
               env={}, cwd=work_dir, check=True)
subprocess.run(['git', 'reset', '--quiet', '--hard', fake_commit],
               env={}, cwd=work_dir, check=True)


if '--fix' in sys.argv:
    print('Running linter to fix errors...')
    p = subprocess.run([linter, '--fix', '--files'] + files)

    # Discover which files were modified
    p = subprocess.run([linter, '--fix', '--files'] + files, cwd=work_dir)
    p = subprocess.run(['git', 'diff', '--name-only', '--', '*.py', '*.pyx'],
                       capture_output=True, check=True, cwd=work_dir)
    files = p.stdout.decode(sys.getfilesystemencoding()).split('\0')
    files = [f for f in files if f]
    if files:
        print('The following files were modified:')
        print()
        print('\n'.join(files))
    else:
        print('No files were modified.\n')

    print('Please remember to `git add` modified files.')
    sys.exit(p.returncode)


p = subprocess.run([linter, '--files'] + files, cwd=work_dir)

if p.returncode != 0:
    print('!! Linting failed; please fix errors, `git add` files, and re-commit.')
    print()
    print('Some errors may be fixable automatically by running:')
    print()
    print('  ./tools/pre-commit-hook.py --fix')

    sys.exit(p.returncode)

p = subprocess.run(unicode_check, cwd=work_dir)

if p.returncode != 0:
    sys.exit(p.returncode)
