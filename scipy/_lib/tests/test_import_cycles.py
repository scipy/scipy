import sys
import subprocess
import pytest

from .test_public_api import PUBLIC_MODULES, PRIVATE_BUT_PRESENT_MODULES

# Regression tests for gh-6793.
# Check that all modules are importable in a new Python process.
# This is not necessarily true if there are import cycles present.

def test_public_modules_importable():
    for module in PUBLIC_MODULES:
        cmd = f'import {module}'
        subprocess.check_call([sys.executable, '-c', cmd])

@pytest.mark.timeout(180)
def test_private_but_present_modules_importable():
    for module in PRIVATE_BUT_PRESENT_MODULES:
        cmd = f'import {module}'
        subprocess.check_call([sys.executable, '-c', cmd])
