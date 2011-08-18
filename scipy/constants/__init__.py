"""
==================================
Constants (:mod:`scipy.constants`)
==================================

.. module:: scipy.constants

Physical and mathematical constants and units.


Mathematical constants
======================

============  =================================================================
``pi``        Pi
``golden``    Golden ratio
============  =================================================================


Physical constants
==================

=============  =================================================================
``c``          speed of light in vacuum
``mu_0``       the magnetic constant :math:`\mu_0`
``epsilon_0``  the electric constant (vacuum permittivity), :math:`\epsilon_0`
``h``          the Planck constant :math:`h`
``hbar``       :math:`\hbar = h/(2\pi)`
``G``          Newtonian constant of gravitation
``g``          standard acceleration of gravity
``e``          elementary charge
``R``          molar gas constant
``alpha``      fine-structure constant
``N_A``        Avogadro constant
``k``          Boltzmann constant
``sigma``      Stefan-Boltzmann constant :math:`\sigma`
``Wien``       Wien displacement law constant
``Rydberg``    Rydberg constant
``m_e``        electron mass
``m_p``        proton mass
``m_n``        neutron mass
=============  =================================================================


Constants database
------------------

In addition to the above variables, :mod:`scipy.constants` also contains the
2010 CODATA recommended values [CODATA2010]_ database containing more physical
constants.

.. autosummary::
   :toctree: generated/

   value      -- Value in physical_constants indexed by key
   unit       -- Unit in physical_constants indexed by key
   precision  -- Relative precision in physical_constants indexed by key
   find       -- Return list of physical_constant keys with a given string
   ConstantWarning -- Constant sought not in newest CODATA data set

.. data:: physical_constants

   Dictionary of physical constants, of the format
   ``physical_constants[name] = (value, unit, uncertainty)``.

Available constants:

======================================================================  ====
%(constant_names)s
======================================================================  ====


Units
=====

SI prefixes
-----------

============  =================================================================
``yotta``     :math:`10^{24}`
``zetta``     :math:`10^{21}`
``exa``       :math:`10^{18}`
``peta``      :math:`10^{15}`
``tera``      :math:`10^{12}`
``giga``      :math:`10^{9}`
``mega``      :math:`10^{6}`
``kilo``      :math:`10^{3}`
``hecto``     :math:`10^{2}`
``deka``      :math:`10^{1}`
``deci``      :math:`10^{-1}`
``centi``     :math:`10^{-2}`
``milli``     :math:`10^{-3}`
``micro``     :math:`10^{-6}`
``nano``      :math:`10^{-9}`
``pico``      :math:`10^{-12}`
``femto``     :math:`10^{-15}`
``atto``      :math:`10^{-18}`
``zepto``     :math:`10^{-21}`
============  =================================================================

Binary prefixes
---------------

============  =================================================================
``kibi``      :math:`2^{10}`
``mebi``      :math:`2^{20}`
``gibi``      :math:`2^{30}`
``tebi``      :math:`2^{40}`
``pebi``      :math:`2^{50}`
``exbi``      :math:`2^{60}`
``zebi``      :math:`2^{70}`
``yobi``      :math:`2^{80}`
============  =================================================================

Weight
------

=================  ============================================================
``gram``           :math:`10^{-3}` kg
``metric_ton``     :math:`10^{3}` kg
``grain``          one grain in kg
``lb``             one pound (avoirdupous) in kg
``oz``             one ounce in kg
``stone``          one stone in kg
``grain``          one grain in kg
``long_ton``       one long ton in kg
``short_ton``      one short ton in kg
``troy_ounce``     one Troy ounce in kg
``troy_pound``     one Troy pound in kg
``carat``          one carat in kg
``m_u``            atomic mass constant (in kg)
=================  ============================================================

Angle
-----

=================  ============================================================
``degree``         degree in radians
``arcmin``         arc minute in radians
``arcsec``         arc second in radians
=================  ============================================================


Time
----

=================  ============================================================
``minute``         one minute in seconds
``hour``           one hour in seconds
``day``            one day in seconds
``week``           one week in seconds
``year``           one year (365 days) in seconds
``Julian_year``    one Julian year (365.25 days) in seconds
=================  ============================================================


Length
------

=================  ============================================================
``inch``           one inch in meters
``foot``           one foot in meters
``yard``           one yard in meters
``mile``           one mile in meters
``mil``            one mil in meters
``pt``             one point in meters
``survey_foot``    one survey foot in meters
``survey_mile``    one survey mile in meters
``nautical_mile``  one nautical mile in meters
``fermi``          one Fermi in meters
``angstrom``       one Angstrom in meters
``micron``         one micron in meters
``au``             one astronomical unit in meters
``light_year``     one light year in meters
``parsec``         one parsec in meters
=================  ============================================================

Pressure
--------

=================  ============================================================
``atm``            standard atmosphere in pascals
``bar``            one bar in pascals
``torr``           one torr (mmHg) in pascals
``psi``            one psi in pascals
=================  ============================================================

Area
----

=================  ============================================================
``hectare``        one hectare in square meters
``acre``           one acre in square meters
=================  ============================================================


Volume
------

===================    ========================================================
``liter``              one liter in cubic meters
``gallon``             one gallon (US) in cubic meters
``gallon_imp``         one gallon (UK) in cubic meters
``fluid_ounce``        one fluid ounce (US) in cubic meters
``fluid_ounce_imp``    one fluid ounce (UK) in cubic meters
``bbl``                one barrel in cubic meters
===================    ========================================================

Speed
-----

=================    ==========================================================
``kmh``              kilometers per hour in meters per second
``mph``              miles per hour in meters per second
``mach``             one Mach (approx., at 15 C, 1 atm) in meters per second
``knot``             one knot in meters per second
=================    ==========================================================


Temperature
-----------

=====================  =======================================================
``zero_Celsius``       zero of Celsius scale in Kelvin
``degree_Fahrenheit``  one Fahrenheit (only differences) in Kelvins
=====================  =======================================================

.. autosummary::
   :toctree: generated/

   C2K
   K2C
   F2C
   C2F
   F2K
   K2F

Energy
------

====================  =======================================================
``eV``                one electron volt in Joules
``calorie``           one calorie (thermochemical) in Joules
``calorie_IT``        one calorie (International Steam Table calorie, 1956) in Joules
``erg``               one erg in Joules
``Btu``               one British thermal unit (International Steam Table) in Joules
``Btu_th``            one British thermal unit (thermochemical) in Joules
``ton_TNT``           one ton of TNT in Joules
====================  =======================================================

Power
-----

====================  =======================================================
``hp``                one horsepower in watts
====================  =======================================================

Force
-----

====================  =======================================================
``dyn``               one dyne in newtons
``lbf``               one pound force in newtons
``kgf``               one kilogram force in newtons
====================  =======================================================

Optics
------

.. autosummary::
   :toctree: generated/

   lambda2nu
   nu2lambda

References
==========

.. [CODATA2010] CODATA Recommended Values of the Fundamental
   Physical Constants 2010.

   http://physics.nist.gov/cuu/Constants/index.html

"""

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
