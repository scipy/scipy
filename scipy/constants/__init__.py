r"""
==================================
Constants (:mod:`scipy.constants`)
==================================

.. currentmodule:: scipy.constants

``scipy.constants.physical_constants[name] = (value, unit, uncertainty)``
contains the 2018 self-consistent values of constants and conversion factors
of physics and chemistry recommended by the Committee on Data of the
International Science Council (CODATA) [JPCRD] [RMP].

They can be searched by ``find`` and are also listed below.
   
.. data:: physical_constants

======================================================================  ====
%(constant_names)s
======================================================================  ====

References
==========

.. [JPCRD] E. Tiesinga, P. J. Mohr, D. B. Newell, and B. N. Taylor, J. Phys.
           Chem. Ref. Data 50, 033105 (2021)
.. [RMP] E. Tiesinga, P. J. Mohr, D. B. Newell, and B. N. Taylor, Rev. Mod.
         Phys. 93, 025010 (2021)

"""
# Modules contributed by BasSw (wegwerp@gmail.com)
from ._codata import *
from ._codata import _obsolete_constants, physical_constants

_constant_names_list = [(_k.lower(), _k, _v)
                        for _k, _v in physical_constants.items()
                        if _k not in _obsolete_constants]
_constant_names = "\n".join(["``%s``%s  %s %s" % (_x[1], " "*(66-len(_x[1])),
                                                  _x[2][0], _x[2][1])
                             for _x in sorted(_constant_names_list)])
if __doc__:
    __doc__ = __doc__ % dict(constant_names=_constant_names)

del _constant_names
del _constant_names_list

__all__ = [s for s in dir() if not s.startswith('_')]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
