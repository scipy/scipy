#! /usr/bin/env python
# Last Change: Mon Aug 20 02:00 PM 2007 J
from __future__ import division, print_function, absolute_import

try:
    from functools import partial
except ImportError:
    from .myfunctools import partial
