#!/bin/bash
#
# Pre-commit linting hook.
#
# Install from root of repository with:
#
#   ln -s tools/pre-commit-hook.sh .git/hooks/pre-commit

# store names of files that were staged
mapfile -d '' changed < <(git diff --cached --name-only -z --diff-filter=ACMR -- '*.py' '*.pyx' '*.pxd' '*.pxi')

# create a temporary copy of what would get committed, without unstaged modifications
# (e.g., only certain changes in a file may have been committed)
GIT_DIR=${GIT_DIR:-.git}
workdir=$GIT_DIR/.pre-commit-workdir
fakecommit=$(git commit-tree -p HEAD $(git write-tree) -m ...)

if ! [[ -d $workdir ]]; then
    git clone -qns "$GIT_DIR" "$workdir" || exit
fi

unset "${!GIT_@}"
cd "$workdir"
git checkout -q $fakecommit

# run linter(s); if any fail, remember that, and set an overall hook failure status
ret=0

# Run lint.py from the scipy source tree
../../tools/lint.py --files "${changed[@]}" || ret=1

if [[ $ret -ne 0 ]]; then
    echo "!! Linting failed; please make fixes, \`git add\` files, and re-commit."
fi
exit $ret
