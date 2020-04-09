#!/usr/bin/env bash

# Remove any existing distribution archives
rm -rf dist
mkdir dist

# Refresh the git hash
git log --pretty=format:'%h' -n 1  > GITHASH

# Generate distribution archives (both source and binary)
pip install --upgrade setuptools wheel
python setup.py sdist
#python setup.py bdist_wheel

# Upload
pip install --upgrade twine
python -m twine upload dist/*
