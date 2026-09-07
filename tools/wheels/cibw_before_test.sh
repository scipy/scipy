set -xe

PROJECT_DIR="${1:-$PWD}"
SCIPY_SRC_DIR="${1:-$PWD}"

# install test dependencies via uv
# the lock file lives in this repo, not in the scipy checkout
PYTHON_EXE="$(python -c 'import sys; print(sys.executable)')"
uv export --project "$PROJECT_DIR" --only-group test-core --frozen | \
    uv pip install --python "$PYTHON_EXE" --no-deps --require-hashes -r -
