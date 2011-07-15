from info import __doc__

# Modules contributed by BasSw (wegwerp@gmail.com)
from codata import *
from constants import *
from codata import _obsolete_constants

_constant_names = [(_k.lower(), _k, _v)
                   for _k, _v in physical_constants.items()
                   if _k not in _obsolete_constants]
_constant_names = "\n".join(["``%s``%s  %s %s" % (_x[1], " "*(66-len(_x[1])),
                                                  _x[2][0], _x[2][1])
                             for _x in sorted(_constant_names)])
__doc__ = __doc__ % dict(constant_names=_constant_names)
del _constant_names

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
