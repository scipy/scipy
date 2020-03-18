r"""
==================================
Constants (:mod:`scipy.constants`)
==================================

.. currentmodule:: scipy.constants


About the SI units
==================
In 2019 seven système international d’unités (SI) base units were defined to
serve as a basis for all physical quantities [CompleteBrochure]. They can be
summarized as followed even though for a detailed background the given
literature is recommended [Intro, ConciseBrochure, CompleteBrochure]:

- The second :math:`\mathrm{s}` as the SI unit of time is derived from the
unperturbed ground-state hyperfine transition frequency of the caesium-133 atom
:math:`\Delta \nu_{\mathrm{Cs}}`.
- The metre :math:`\mathrm{m}` as the SI unit of length is derived from the to
the speed of light in vaccum :math:`c`.
- The kilogram :math:`\mathrm{kg}` as the SI unit of mass is derived from the
Planck constant :math:`h`.
- The ampere :math:`\mathrm{A}` as the SI unit of electric current is derived
from the elementary charge :math:`e`.
- The kelvin :math:`\mathrm{K}` as the SI unit of thermodynamic temperature is
derived from the Boltzmann constant :math:`k`.
- The mole :math:`\mathrm{mol}` as the SI unit of amount of substance is
derived from Avogadro constant :math:`N_{\mathrm{A}}`.
- The candela :math:`\mathrm{cd}` as the SI unit of luminous intensity in a
given direction is derived from the luminous efficacy of monochromatic
radiation of frequency :math:`K_{\mathrm{cd}}`.

Defining these SI base units as constants has the advantage will remain for all
future as they are defined now and all physical quantities can expressed by a
combination of these base units. Also all "traditional" units such as newton or
pascal can be expressed in SI units. Therefore the SI brochure list 22 derived
units with special names in the SI and serveral non-SI units. Of course all of
before mentioned units can be used with the SI prefixes.


About the notation of physical units
====================================
A physical quantity is written as the product of a number and a unit such as in
:math:`T = 293 \ \mathrm{K}` where :math:`293` is the number and
:math:`\mathrm{K}` is the unit. The quantity :math:`T` is written in italic
(slanted) symbols and the unit :math:`\mathrm{K}` is written in roman (upright)
symbols. As quantities are written as a product of a number and a unit they are
treated by the common rules of algebra meaning that
:math:`T = 293 \ \mathrm{K}` equals :math:`T/\mathrm{K} = 293`. The latter form
is also the recommended form for tables or figures where :math:`T/\mathrm{K}`
is put in the table head or on the axis label and the numerical values are put
in the table rows or the axis ticks.


List of all constants
=====================
The physical constants provided by SciPy are given as
``physical_constants[name] = (value, unit, uncertainty)``
where
   value      -- Value in physical_constants indexed by key
   unit       -- Unit in physical_constants indexed by key
   precision  -- Relative precision in physical_constants indexed by key
   find       -- Return list of physical_constant keys with a given string
   ConstantWarning -- Constant sought not in newest CODATA data set

The complete list of available constants is given in the following:
%(constant_names)s


References
==========
.. [Intro] Introduction to SI units, https://www.bipm.org/en/measurement-units/
.. [CompleteBrochure] Le Système international d’unités (SI), 9e édition 2019,
   v1.06, https://www.bipm.org/utils/common/pdf/si-brochure/SI-Brochure-9.pdf
.. [ConciseBrochure] A concise summary of theInternational System of Units, SI,
   https://www.bipm.org/utils/common/pdf/si-brochure/SI-Brochure-9-concise-EN.pdf
.. [Wikipedia] Wikipedia "International System of Units",
   https://en.wikipedia.org/wiki/International_System_of_Units
.. [CODATA2018] Fundamental Physical Constants --- Complete Listing, 2018
   CODATA adjustment 2018, https://physics.nist.gov/cuu/Constants/

# Modules contributed by BasSw (wegwerp@gmail.com)
from .codata import *
from .constants import *
from .codata import _obsolete_constants

_constant_names = [(_k.lower(), _k, _v)
                   for _k, _v in physical_constants.items()
                   if _k not in _obsolete_constants]
_constant_names = "\n".join(["``%s``%s  %s %s" % (_x[1], " "*(66-len(_x[1])),
                                                  _x[2][0], _x[2][1])
                             for _x in sorted(_constant_names)])
if __doc__:
    __doc__ = __doc__ % dict(constant_names=_constant_names)

del _constant_names

__all__ = [s for s in dir() if not s.startswith('_')]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
