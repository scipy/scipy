# Import correct md5, irrespective of the Python version
#
# `hashlib` was introduced in 2.5, deprecating `md5`

try:
    from hashlib import *
except ImportError:
    from md5 import *

new = md5
