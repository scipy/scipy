#
# stats - Statistical Functions
#

from info_stats import __doc__

from stats import *
from distributions import *
from rv import *
from morestats import *
from kde import gaussian_kde

try:  # use R functions if installed.
    import rpy
    from rfuncs import *
except ImportError:
    pass

