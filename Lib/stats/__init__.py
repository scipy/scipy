#
# stats - Statistical Functions
#

from info_stats import __doc__

import pstats
from stats import *
from distributions import *
from rv import *
from morestats import *

try:  # use R functions if installed.
    import rpy
    from rfuncs import *
except ImportError:
    pass

