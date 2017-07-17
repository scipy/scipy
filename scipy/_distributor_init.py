""" Distributor init file

Distributors: you can add custom code here to support particular distributions
of scipy.

For example, this is a good place to put any checks for hardware requirements.

The scipy standard source distribution will not put code in this file, so you
can safely replace this file with your own version.
"""

import os
import sys

# Patch the DLL search path on windows so that extensions can load any DLLs that
# they depend on (specifically gfortran DLLs) on windows.

if sys.platform == 'win32':
    os.environ["PATH"] += os.pathsep + os.path.join(os.path.dirname(__file__), '_lib')
