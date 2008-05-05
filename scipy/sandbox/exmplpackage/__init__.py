#
# exmplpackage - Example of a Scipy module, see DEVELOPERS.txt
#

# Get documentation string:
from info_exmplpackage import __doc__

# Import testing rig, allowing scipy.examplpackage.test()
from scipy.testing.pkgtester import Tester
test = Tester().test

# Import symbols from sub-module:
from foo import *

# Import sub-package:
import yyy

# Import extension module
import spam
