import math
import subprocess
import sys

import pytest

from .test_public_api import PUBLIC_MODULES

# Regression tests for gh-6793.
# Check that all modules are importable in a new Python process.
# This is not necessarily true if there are import cycles present.

@pytest.mark.fail_slow(40)
@pytest.mark.slow
def test_public_modules_importable():
    # Split into batches to limit peak resource usage (memory, file handles)
    # on resource-constrained systems (e.g., RISC-V - see gh-24163).
    # A regular for-loop over all modules is too slow (~4x slower).
    n_batches = 6
    batch_size = math.ceil(len(PUBLIC_MODULES) / n_batches)

    for i in range(n_batches):
        batch = PUBLIC_MODULES[i * batch_size:(i + 1) * batch_size]
        pids = [subprocess.Popen([sys.executable, '-c', f'import {module}'])
                for module in batch]
        for j, pid in enumerate(pids):
            assert pid.wait() == 0, f'Failed to import {batch[j]}'
