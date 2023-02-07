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

subprocess.run(['git', 'checkout', '-q', fake_commit],
               env={}, cwd=work_dir, check=True)

# Run lint.py from the scipy source tree
linters = [
    '../../tools/lint.py',
    'tools/lint.py',
    'lint.py'  # in case pre-commit hook is run from tools dir
]

linter = [f for f in linters if os.path.exists(f)][0]

p = subprocess.run([linter, '--fix', '--files'] + files,
                   capture_output=True, text=True)
print(p.stdout.strip(), end='')
print(p.stderr.strip(), end='')

if 'fixed, 0 remaining' in p.stdout:
    print("\n\nAll errors have been fixed; please `git add` and re-commit.")
    sys.exit(1)

if p.returncode != 0:
    print("\n\n!! Linting failed; please make fixes, `git add` files, and re-commit.")
    sys.exit(p.returncode)
