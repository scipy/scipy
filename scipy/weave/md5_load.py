# Import correct md5, irrespective of the Python version
#
# `hashlib` was introduced in 2.5, deprecating `md5`
from __future__ import absolute_import, print_function

try:
    from hashlib import *
except ImportError:
    from md5 import *

new = md5
