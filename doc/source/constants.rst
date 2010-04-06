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
==================

In addition to the above variables containing physical constants,
:mod:`scipy.constants` also contains a database of additional physical
constants.

.. autosummary::
   :toctree: generated/

   value
   unit
   precision
   find

.. data:: physical_constants

   Dictionary of physical constants, of the format
   ``physical_constants[name] = (value, unit, uncertainty)``.

   Available constants:

   ======================================================================  ====
   ``alpha particle mass``
   ``alpha particle mass energy equivalent``
   ``alpha particle mass energy equivalent in MeV``
   ``alpha particle mass in u``
   ``alpha particle molar mass``
   ``alpha particle-electron mass ratio``
   ``alpha particle-proton mass ratio``
   ``Angstrom star``
   ``atomic mass constant``
   ``atomic mass constant energy equivalent``
   ``atomic mass constant energy equivalent in MeV``
   ``atomic mass unit-electron volt relationship``
   ``atomic mass unit-hartree relationship``
   ``atomic mass unit-hertz relationship``
   ``atomic mass unit-inverse meter relationship``
   ``atomic mass unit-joule relationship``
   ``atomic mass unit-kelvin relationship``
   ``atomic mass unit-kilogram relationship``
   ``atomic unit of 1st hyperpolarizablity``
   ``atomic unit of 2nd hyperpolarizablity``
   ``atomic unit of action``
   ``atomic unit of charge``
   ``atomic unit of charge density``
   ``atomic unit of current``
   ``atomic unit of electric dipole moment``
   ``atomic unit of electric field``
   ``atomic unit of electric field gradient``
   ``atomic unit of electric polarizablity``
   ``atomic unit of electric potential``
   ``atomic unit of electric quadrupole moment``
   ``atomic unit of energy``
   ``atomic unit of force``
   ``atomic unit of length``
   ``atomic unit of magnetic dipole moment``
   ``atomic unit of magnetic flux density``
   ``atomic unit of magnetizability``
   ``atomic unit of mass``
   ``atomic unit of momentum``
   ``atomic unit of permittivity``
   ``atomic unit of time``
   ``atomic unit of velocity``
   ``Avogadro constant``
   ``Bohr magneton``
   ``Bohr magneton in eV/T``
   ``Bohr magneton in Hz/T``
   ``Bohr magneton in inverse meters per tesla``
   ``Bohr magneton in K/T``
   ``Bohr radius``
   ``Boltzmann constant``
   ``Boltzmann constant in eV/K``
   ``Boltzmann constant in Hz/K``
   ``Boltzmann constant in inverse meters per kelvin``
   ``characteristic impedance of vacuum``
   ``classical electron radius``
   ``Compton wavelength``
   ``Compton wavelength over 2 pi``
   ``conductance quantum``
   ``conventional value of Josephson constant``
   ``conventional value of von Klitzing constant``
   ``Cu x unit``
   ``deuteron magnetic moment``
   ``deuteron magnetic moment to Bohr magneton ratio``
   ``deuteron magnetic moment to nuclear magneton ratio``
   ``deuteron mass``
   ``deuteron mass energy equivalent``
   ``deuteron mass energy equivalent in MeV``
   ``deuteron mass in u``
   ``deuteron molar mass``
   ``deuteron rms charge radius``
   ``deuteron-electron magnetic moment ratio``
   ``deuteron-electron mass ratio``
   ``deuteron-neutron magnetic moment ratio``
   ``deuteron-proton magnetic moment ratio``
   ``deuteron-proton mass ratio``
   ``electric constant``
   ``electron charge to mass quotient``
   ``electron g factor``
   ``electron gyromagnetic ratio``
   ``electron gyromagnetic ratio over 2 pi``
   ``electron magnetic moment``
   ``electron magnetic moment anomaly``
   ``electron magnetic moment to Bohr magneton ratio``
   ``electron magnetic moment to nuclear magneton ratio``
   ``electron mass``
   ``electron mass energy equivalent``
   ``electron mass energy equivalent in MeV``
   ``electron mass in u``
   ``electron molar mass``
   ``electron to alpha particle mass ratio``
   ``electron to shielded helion magnetic moment ratio``
   ``electron to shielded proton magnetic moment ratio``
   ``electron volt``
   ``electron volt-atomic mass unit relationship``
   ``electron volt-hartree relationship``
   ``electron volt-hertz relationship``
   ``electron volt-inverse meter relationship``
   ``electron volt-joule relationship``
   ``electron volt-kelvin relationship``
   ``electron volt-kilogram relationship``
   ``electron-deuteron magnetic moment ratio``
   ``electron-deuteron mass ratio``
   ``electron-muon magnetic moment ratio``
   ``electron-muon mass ratio``
   ``electron-neutron magnetic moment ratio``
   ``electron-neutron mass ratio``
   ``electron-proton magnetic moment ratio``
   ``electron-proton mass ratio``
   ``electron-tau mass ratio``
   ``elementary charge``
   ``elementary charge over h``
   ``Faraday constant``
   ``Faraday constant for conventional electric current``
   ``Fermi coupling constant``
   ``fine-structure constant``
   ``first radiation constant``
   ``first radiation constant for spectral radiance``
   ``Hartree energy``
   ``Hartree energy in eV``
   ``hartree-atomic mass unit relationship``
   ``hartree-electron volt relationship``
   ``hartree-hertz relationship``
   ``hartree-inverse meter relationship``
   ``hartree-joule relationship``
   ``hartree-kelvin relationship``
   ``hartree-kilogram relationship``
   ``helion mass``
   ``helion mass energy equivalent``
   ``helion mass energy equivalent in MeV``
   ``helion mass in u``
   ``helion molar mass``
   ``helion-electron mass ratio``
   ``helion-proton mass ratio``
   ``hertz-atomic mass unit relationship``
   ``hertz-electron volt relationship``
   ``hertz-hartree relationship``
   ``hertz-inverse meter relationship``
   ``hertz-joule relationship``
   ``hertz-kelvin relationship``
   ``hertz-kilogram relationship``
   ``inverse fine-structure constant``
   ``inverse meter-atomic mass unit relationship``
   ``inverse meter-electron volt relationship``
   ``inverse meter-hartree relationship``
   ``inverse meter-hertz relationship``
   ``inverse meter-joule relationship``
   ``inverse meter-kelvin relationship``
   ``inverse meter-kilogram relationship``
   ``inverse of conductance quantum``
   ``Josephson constant``
   ``joule-atomic mass unit relationship``
   ``joule-electron volt relationship``
   ``joule-hartree relationship``
   ``joule-hertz relationship``
   ``joule-inverse meter relationship``
   ``joule-kelvin relationship``
   ``joule-kilogram relationship``
   ``kelvin-atomic mass unit relationship``
   ``kelvin-electron volt relationship``
   ``kelvin-hartree relationship``
   ``kelvin-hertz relationship``
   ``kelvin-inverse meter relationship``
   ``kelvin-joule relationship``
   ``kelvin-kilogram relationship``
   ``kilogram-atomic mass unit relationship``
   ``kilogram-electron volt relationship``
   ``kilogram-hartree relationship``
   ``kilogram-hertz relationship``
   ``kilogram-inverse meter relationship``
   ``kilogram-joule relationship``
   ``kilogram-kelvin relationship``
   ``lattice parameter of silicon``
   ``Loschmidt constant (273.15 K, 101.325 kPa)``
   ``magnetic constant``
   ``magnetic flux quantum``
   ``Mo x unit``
   ``molar gas constant``
   ``molar mass constant``
   ``molar mass of carbon-12``
   ``molar Planck constant``
   ``molar Planck constant times c``
   ``molar volume of ideal gas (273.15 K, 100 kPa)``
   ``molar volume of ideal gas (273.15 K, 101.325 kPa)``
   ``molar volume of silicon``
   ``muon Compton wavelength``
   ``muon Compton wavelength over 2 pi``
   ``muon g factor``
   ``muon magnetic moment``
   ``muon magnetic moment anomaly``
   ``muon magnetic moment to Bohr magneton ratio``
   ``muon magnetic moment to nuclear magneton ratio``
   ``muon mass``
   ``muon mass energy equivalent``
   ``muon mass energy equivalent in MeV``
   ``muon mass in u``
   ``muon molar mass``
   ``muon-electron mass ratio``
   ``muon-neutron mass ratio``
   ``muon-proton magnetic moment ratio``
   ``muon-proton mass ratio``
   ``muon-tau mass ratio``
   ``natural unit of action``
   ``natural unit of action in eV s``
   ``natural unit of energy``
   ``natural unit of energy in MeV``
   ``natural unit of length``
   ``natural unit of mass``
   ``natural unit of momentum``
   ``natural unit of momentum in MeV/c``
   ``natural unit of time``
   ``natural unit of velocity``
   ``neutron Compton wavelength``
   ``neutron Compton wavelength over 2 pi``
   ``neutron g factor``
   ``neutron gyromagnetic ratio``
   ``neutron gyromagnetic ratio over 2 pi``
   ``neutron magnetic moment``
   ``neutron magnetic moment to Bohr magneton ratio``
   ``neutron magnetic moment to nuclear magneton ratio``
   ``neutron mass``
   ``neutron mass energy equivalent``
   ``neutron mass energy equivalent in MeV``
   ``neutron mass in u``
   ``neutron molar mass``
   ``neutron to shielded proton magnetic moment ratio``
   ``neutron-electron magnetic moment ratio``
   ``neutron-electron mass ratio``
   ``neutron-muon mass ratio``
   ``neutron-proton magnetic moment ratio``
   ``neutron-proton mass ratio``
   ``neutron-tau mass ratio``
   ``Newtonian constant of gravitation``
   ``Newtonian constant of gravitation over h-bar c``
   ``nuclear magneton``
   ``nuclear magneton in eV/T``
   ``nuclear magneton in inverse meters per tesla``
   ``nuclear magneton in K/T``
   ``nuclear magneton in MHz/T``
   ``Planck constant``
   ``Planck constant in eV s``
   ``Planck constant over 2 pi``
   ``Planck constant over 2 pi in eV s``
   ``Planck constant over 2 pi times c in MeV fm``
   ``Planck length``
   ``Planck mass``
   ``Planck temperature``
   ``Planck time``
   ``proton charge to mass quotient``
   ``proton Compton wavelength``
   ``proton Compton wavelength over 2 pi``
   ``proton g factor``
   ``proton gyromagnetic ratio``
   ``proton gyromagnetic ratio over 2 pi``
   ``proton magnetic moment``
   ``proton magnetic moment to Bohr magneton ratio``
   ``proton magnetic moment to nuclear magneton ratio``
   ``proton magnetic shielding correction``
   ``proton mass``
   ``proton mass energy equivalent``
   ``proton mass energy equivalent in MeV``
   ``proton mass in u``
   ``proton molar mass``
   ``proton rms charge radius``
   ``proton-electron mass ratio``
   ``proton-muon mass ratio``
   ``proton-neutron magnetic moment ratio``
   ``proton-neutron mass ratio``
   ``proton-tau mass ratio``
   ``quantum of circulation``
   ``quantum of circulation times 2``
   ``Rydberg constant``
   ``Rydberg constant times c in Hz``
   ``Rydberg constant times hc in eV``
   ``Rydberg constant times hc in J``
   ``Sackur-Tetrode constant (1 K, 100 kPa)``
   ``Sackur-Tetrode constant (1 K, 101.325 kPa)``
   ``second radiation constant``
   ``shielded helion gyromagnetic ratio``
   ``shielded helion gyromagnetic ratio over 2 pi``
   ``shielded helion magnetic moment``
   ``shielded helion magnetic moment to Bohr magneton ratio``
   ``shielded helion magnetic moment to nuclear magneton ratio``
   ``shielded helion to proton magnetic moment ratio``
   ``shielded helion to shielded proton magnetic moment ratio``
   ``shielded proton gyromagnetic ratio``
   ``shielded proton gyromagnetic ratio over 2 pi``
   ``shielded proton magnetic moment``
   ``shielded proton magnetic moment to Bohr magneton ratio``
   ``shielded proton magnetic moment to nuclear magneton ratio``
   ``speed of light in vacuum``
   ``standard acceleration of gravity``
   ``standard atmosphere``
   ``Stefan-Boltzmann constant``
   ``tau Compton wavelength``
   ``tau Compton wavelength over 2 pi``
   ``tau mass``
   ``tau mass energy equivalent``
   ``tau mass energy equivalent in MeV``
   ``tau mass in u``
   ``tau molar mass``
   ``tau-electron mass ratio``
   ``tau-muon mass ratio``
   ``tau-neutron mass ratio``
   ``tau-proton mass ratio``
   ``Thomson cross section``
   ``unified atomic mass unit``
   ``von Klitzing constant``
   ``weak mixing angle``
   ``Wien displacement law constant``
   ``{220} lattice spacing of silicon``
   ======================================================================  ====


Unit prefixes
=============

SI
--

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


Binary
------

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

Units
=====

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
``angstrom``       one Ångström in meters
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
``mach``             one Mach (approx., at 15 °C, 1 atm) in meters per second
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
