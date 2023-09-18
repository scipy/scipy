# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys

# This __version__ assignment is parsed by setup.py; keep it in this form.
# Development versions end with ".dev" (suffix is added below).
__version__ = '0.4.2'
__release__ = not __version__.endswith(".dev")

try:
    from ._version import __githash__, __suffix__
except ImportError:
    __githash__ = None
    __suffix__ = "0"
if not __release__:
    __version__ += __suffix__
del __suffix__


if sys.version_info >= (3, 3):
    # OS X framework builds of Python 3.3 can not call other 3.3
    # virtualenvs as a subprocess because `__PYENV_LAUNCHER__` is
    # inherited.
    if os.environ.get('__PYVENV_LAUNCHER__'):
        os.unsetenv('__PYVENV_LAUNCHER__')

from . import plugin_manager
