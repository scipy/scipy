#
# plt - Plotting with wxPython
#

from info_plt import __doc__

from lena import lena

from plot_utility import *
from interface import *

# XXX: why we need this?:
import plot_utility
import interface

from scipy_base import display_test
if not display_test.have_x11() or display_test.try_XOpenDisplay():
    from plot_objects import *
    from wxplt import *
    import wxplt
    import plot_objects


