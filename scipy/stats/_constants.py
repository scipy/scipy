"""
Statistics-related constants.

"""
from __future__ import division, print_function, absolute_import

import numpy as np


# The smallest representable positive number such that 1.0 + _EPS != 1.0.
_EPS = np.finfo(float).eps

# The largest [in magnitude] usable floating value.
_XMAX = np.finfo(float).machar.xmax

# The smallest [in magnitude] usable floating value.
_XMIN = np.finfo(float).machar.xmin

# -special.psi(1)
_EULER = 0.577215664901532860606512090082402431042

# special.zeta(3, 1)  Apery's constant
_ZETA3 = 1.202056903159594285399738161511449990765
