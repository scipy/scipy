from info import __doc__

# Modules contributed by BasSw (wegwerp@gmail.com)
from codata import *
from constants import *

_constant_names = [(k.lower(), k, v) for k, v in physical_constants.items()]
_constant_names = "\n".join(["``%s``%s  %s %s" % (x[1], " "*(66-len(x[1])),
                                                  x[2][0], x[2][1])
                             for x in sorted(_constant_names)])
__doc__ = __doc__ % dict(constant_names=_constant_names)
del _constant_names

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
