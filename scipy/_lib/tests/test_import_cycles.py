import sys
import subprocess

from .test_public_api import PUBLIC_MODULES

# Regression tests for gh-6793.
# Check that all modules are importable in a new Python process.
# This is not necessarily true if there are import cycles present.

def test_public_modules_importable():
    for module in PUBLIC_MODULES:
        cmd = f'import {module}'
        subprocess.check_call([sys.executable, '-c', cmd])
