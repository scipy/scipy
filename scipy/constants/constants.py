"""
Collection of physical constants and conversion factors.

Most constants are in SI units, so you can do
print '10 mile per minute is', 10*mile/minute, 'm/s or', 10*mile/(minute*knot), 'knots'

The list is not meant to be comprehensive, but just a convenient list for everyday use.
"""
from __future__ import division, print_function, absolute_import

"""
BasSw 2006
physical constants: imported from CODATA
unit conversion: see e.g. NIST special publication 811
Use at own risk: double-check values before calculating your Mars orbit-insertion burn.
Some constants exist in a few variants, which are marked with suffixes.
The ones without any suffix should be the most common one.
"""

import math as _math
from .codata import value as _cd
import numpy as _np

# mathematical constants
pi = _math.pi
golden = golden_ratio = (1 + _math.sqrt(5)) / 2

# SI prefixes
yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1e3
hecto = 1e2
deka = 1e1
deci = 1e-1
centi = 1e-2
milli = 1e-3
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21

# binary prefixes
kibi = 2**10
mebi = 2**20
gibi = 2**30
tebi = 2**40
pebi = 2**50
exbi = 2**60
zebi = 2**70
yobi = 2**80

# physical constants
c = speed_of_light = _cd('speed of light in vacuum')
mu_0 = 4e-7*pi
epsilon_0 = 1 / (mu_0*c*c)
h = Planck = _cd('Planck constant')
hbar = h / (2 * pi)
G = gravitational_constant = _cd('Newtonian constant of gravitation')
g = _cd('standard acceleration of gravity')
e = elementary_charge = _cd('elementary charge')
R = gas_constant = _cd('molar gas constant')
alpha = fine_structure = _cd('fine-structure constant')
N_A = Avogadro = _cd('Avogadro constant')
k = Boltzmann = _cd('Boltzmann constant')
sigma = Stefan_Boltzmann = _cd('Stefan-Boltzmann constant')
Wien = _cd('Wien wavelength displacement law constant')
Rydberg = _cd('Rydberg constant')

# weight in kg
gram = 1e-3
metric_ton = 1e3
grain = 64.79891e-6
lb = pound = 7000 * grain  # avoirdupois
oz = ounce = pound / 16
stone = 14 * pound
long_ton = 2240 * pound
short_ton = 2000 * pound

troy_ounce = 480 * grain  # only for metals / gems
troy_pound = 12 * troy_ounce
carat = 200e-6

m_e = electron_mass = _cd('electron mass')
m_p = proton_mass = _cd('proton mass')
m_n = neutron_mass = _cd('neutron mass')
m_u = u = atomic_mass = _cd('atomic mass constant')

# angle in rad
degree = pi / 180
arcmin = arcminute = degree / 60
arcsec = arcsecond = arcmin / 60

# time in second
minute = 60.0
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day
Julian_year = 365.25 * day

# length in meter
inch = 0.0254
foot = 12 * inch
yard = 3 * foot
mile = 1760 * yard
mil = inch / 1000
pt = point = inch / 72  # typography
survey_foot = 1200.0 / 3937
survey_mile = 5280 * survey_foot
nautical_mile = 1852.0
fermi = 1e-15
angstrom = 1e-10
micron = 1e-6
au = astronomical_unit = 149597870691.0
light_year = Julian_year * c
parsec = au / arcsec

# pressure in pascal
atm = atmosphere = _cd('standard atmosphere')
bar = 1e5
torr = mmHg = atm / 760
psi = pound * g / (inch * inch)

# area in meter**2
hectare = 1e4
acre = 43560 * foot**2

# volume in meter**3
litre = liter = 1e-3
gallon = gallon_US = 231 * inch**3  # US
# pint = gallon_US / 8
fluid_ounce = fluid_ounce_US = gallon_US / 128
bbl = barrel = 42 * gallon_US  # for oil

gallon_imp = 4.54609e-3  # UK
fluid_ounce_imp = gallon_imp / 160

# speed in meter per second
kmh = 1e3 / hour
mph = mile / hour
mach = speed_of_sound = 340.5  # approx value at 15 degrees in 1 atm. is this a common value?
knot = nautical_mile / hour

# temperature in kelvin
zero_Celsius = 273.15
degree_Fahrenheit = 1/1.8  # only for differences

# energy in joule
eV = electron_volt = elementary_charge  # * 1 Volt
calorie = calorie_th = 4.184
calorie_IT = 4.1868
erg = 1e-7
Btu_th = pound * degree_Fahrenheit * calorie_th / gram
Btu = Btu_IT = pound * degree_Fahrenheit * calorie_IT / gram
ton_TNT = 1e9 * calorie_th
# Wh = watt_hour

# power in watt
hp = horsepower = 550 * foot * pound * g

# force in newton
dyn = dyne = 1e-5
lbf = pound_force = pound * g
kgf = kilogram_force = g  # * 1 kg

# functions for conversions that are not linear


def C2K(C):
    """
    Convert Celsius to Kelvin

    Parameters
    ----------
    C : array_like
        Celsius temperature(s) to be converted.

    Returns
    -------
    K : float or array of floats
        Equivalent Kelvin temperature(s).

    Notes
    -----
    Computes ``K = C + zero_Celsius`` where `zero_Celsius` = 273.15, i.e.,
    (the absolute value of) temperature "absolute zero" as measured in Celsius.

    Examples
    --------
    >>> from scipy.constants import C2K
    >>> C2K(np.array([-40, 40.0]))
    array([ 233.15,  313.15])

    """
    return _np.asanyarray(C) + zero_Celsius


def K2C(K):
    """
    Convert Kelvin to Celsius

    Parameters
    ----------
    K : array_like
        Kelvin temperature(s) to be converted.

    Returns
    -------
    C : float or array of floats
        Equivalent Celsius temperature(s).

    Notes
    -----
    Computes ``C = K - zero_Celsius`` where `zero_Celsius` = 273.15, i.e.,
    (the absolute value of) temperature "absolute zero" as measured in Celsius.

    Examples
    --------
    >>> from scipy.constants import K2C
    >>> K2C(np.array([233.15, 313.15]))
    array([-40.,  40.])

    """
    return _np.asanyarray(K) - zero_Celsius


def F2C(F):
    """
    Convert Fahrenheit to Celsius

    Parameters
    ----------
    F : array_like
        Fahrenheit temperature(s) to be converted.

    Returns
    -------
    C : float or array of floats
        Equivalent Celsius temperature(s).

    Notes
    -----
    Computes ``C = (F - 32) / 1.8``.

    Examples
    --------
    >>> from scipy.constants import F2C
    >>> F2C(np.array([-40, 40.0]))
    array([-40.        ,   4.44444444])

    """
    return (_np.asanyarray(F) - 32) / 1.8


def C2F(C):
    """
    Convert Celsius to Fahrenheit

    Parameters
    ----------
    C : array_like
        Celsius temperature(s) to be converted.

    Returns
    -------
    F : float or array of floats
        Equivalent Fahrenheit temperature(s).

    Notes
    -----
    Computes ``F = 1.8 * C + 32``.

    Examples
    --------
    >>> from scipy.constants import C2F
    >>> C2F(np.array([-40, 40.0]))
    array([ -40.,  104.])

    """
    return 1.8 * _np.asanyarray(C) + 32


def F2K(F):
    """
    Convert Fahrenheit to Kelvin

    Parameters
    ----------
    F : array_like
        Fahrenheit temperature(s) to be converted.

    Returns
    -------
    K : float or array of floats
        Equivalent Kelvin temperature(s).

    Notes
    -----
    Computes ``K = (F - 32)/1.8 + zero_Celsius`` where `zero_Celsius` =
    273.15, i.e., (the absolute value of) temperature "absolute zero" as
    measured in Celsius.

    Examples
    --------
    >>> from scipy.constants import F2K
    >>> F2K(np.array([-40, 104]))
    array([ 233.15,  313.15])

    """
    return C2K(F2C(_np.asanyarray(F)))


def K2F(K):
    """
    Convert Kelvin to Fahrenheit

    Parameters
    ----------
    K : array_like
        Kelvin temperature(s) to be converted.

    Returns
    -------
    F : float or array of floats
        Equivalent Fahrenheit temperature(s).

    Notes
    -----
    Computes ``F = 1.8 * (K - zero_Celsius) + 32`` where `zero_Celsius` =
    273.15, i.e., (the absolute value of) temperature "absolute zero" as
    measured in Celsius.

    Examples
    --------
    >>> from scipy.constants import K2F
    >>> K2F(np.array([233.15,  313.15]))
    array([ -40.,  104.])

    """
    return C2F(K2C(_np.asanyarray(K)))

# optics


def lambda2nu(lambda_):
    """
    Convert wavelength to optical frequency

    Parameters
    ----------
    lambda_ : array_like
        Wavelength(s) to be converted.

    Returns
    -------
    nu : float or array of floats
        Equivalent optical frequency.

    Notes
    -----
    Computes ``nu = c / lambda`` where c = 299792458.0, i.e., the
    (vacuum) speed of light in meters/second.

    Examples
    --------
    >>> from scipy.constants import lambda2nu, speed_of_light
    >>> lambda2nu(np.array((1, speed_of_light)))
    array([  2.99792458e+08,   1.00000000e+00])

    """
    return _np.asanyarray(c) / lambda_


def nu2lambda(nu):
    """
    Convert optical frequency to wavelength.

    Parameters
    ----------
    nu : array_like
        Optical frequency to be converted.

    Returns
    -------
    lambda : float or array of floats
        Equivalent wavelength(s).

    Notes
    -----
    Computes ``lambda = c / nu`` where c = 299792458.0, i.e., the
    (vacuum) speed of light in meters/second.

    Examples
    --------
    >>> from scipy.constants import nu2lambda, speed_of_light
    >>> nu2lambda(np.array((1, speed_of_light)))
    array([  2.99792458e+08,   1.00000000e+00])

    """
    return c / _np.asanyarray(nu)

# length

def ft_to_meters(ft):
    """ 
    Convert foot to meters

    Parameters
    ----------
    ft: array_like
        Foot length to be converted.

    Returns
    -------
    meters: float or array of floats
        Equivalent length in meters.

    Notes 
    -----
    Computes ``meters = ft * 0.3048 ` where 0.3048 meter = 1ft 
    
    """
    return _np.asanyarray(ft)*0.3048

def meters_to_ft(meters):
    """ 
    Convert meters to foot

    Parameters
    ----------
    meters: array_like
        meters length to be converted.

    Returns
    -------
    ft: float or array of floats
        Equivalent length in foot.

    Notes 
    -----
    Computes ``ft = meters/ 0.3048 ` where 0.3048 meter = 1ft 
    
    """
    return _np.asanyarray(meters)/0.3048

def inches_to_ft(in):
    """ 
    Convert foot to inches

    Parameters
    ----------
    ft: array_like
        Foot length to be converted.

    Returns
    -------
    inches: float or array of floats
        Equivalent length in inches.

    Notes 
    -----
    Computes ``in = ft * 12 ` where 12 inches = 1ft 
    
    """
    return (_npasanyarray(ft))*12

def ft_to_inches(ft):
    """ 
    Convert inches to foot

    Parameters
    ----------
    in: array_like
        Inches length to be converted.

    Returns
    -------
    foot: float or array of floats
        Equivalent length in foot.

    Notes 
    -----
    Computes ``ft = in / 12 ` where 12 inches = 1ft 
    
    """
    return (_npasanyarray(in))/12


def in_to_meters(in):
    """ 
    Convert inches to meters

    Parameters
    ----------
    in: array_like
        Inches length to be converted.

    Returns
    -------
    meters: float or array of floats
        Equivalent length in meters.

    Notes 
    -----
    Computes ``meters = 0.0254 * in `` where 1 meter = 0.0254 in.
    
    """
    return 0.0254*(_npasanyarray(in))

def meters_to_in(meters):
    """ 
    Convert meters to inches

    Parameters
    ----------
    meters: array_like
        Meters length to be converted.

    Returns
    -------
    inches: float or array of floats
        Equivalent length in inches.

    Notes 
    -----
    Computes ``inches = meters / 0.0254  `` where 1 meter = 0.0254 in.
    
    """
    return (_npasanyarray(meters))/0.0254

# mass

def lb_to_kg(lb):
    """ 
    Convert foot to meters

    Parameters
    ----------
    lb: array_like
        Pounds to be converted.

    Returns
    -------
    kilogram: float or array of floats
        Equivalent mass in kilogram(s).

    Notes 
    -----
    Computes ``kg = 0.453592 * lb`` where 1 kg = 0.453592 lb.
    
    """
    return 0.453592*_npasanyarray(lb)

def kg_to_lb(kg):
    """ 
    Convert kilograms to pounds 

    Parameters
    ----------
    kg: array_like
        Kilograms to be converted.

    Returns
    -------
    pounds: float or array of floats
        Equivalent mass in pounds.

    Notes 
    -----
    Computes ``lb = kg/0.453592 `` where 1 kg = 0.453592 lb.
    
    """
    return _npasanyarray(kg)/0.453592

def oz_to_g(oz):
    """ 
    Convert ounces to grams

    Parameters
    ----------
    oz: array_like
        Ounces to be converted.

    Returns
    -------
    grams: float or array of floats
        Equivalent masss in grams.

    Notes 
    -----
    Computes ``g = 28.3495 * oz `` where 1 gram = 28.3495 oz.
    
    """
    return 28.3495*_npasanyarray(oz)

def g_to_oz(g):
    """ 
    Convert grams to ounces

    Parameters
    ----------
    g: array_like
        Grams to be converted.

    Returns
    -------
    Ounces: float or array of floats
        Equivalent masss in ounces.

    Notes 
    -----
    Computes ``oz = g / 28.3495 `` where 1 gram = 28.3495 oz.
    
    """
    return _npasanyarray(g) / 28.3495 

# time

def sec_to_min(sec):
    """
    Convert seconds to minutes

    Parameters
    ----------
    sec: array_like
        Seconds to be converted.

    Returns
    -------
    minutes: float or array of floats
        Equivalent time in minutes.

    Notes
    -----
    Computes `` min = sec/60`` where 1 minute = 60 seconds.

    """
    return _npasanyarray(sec)/60

def min_to_sec(min):
    """
    Convert minutes to seconds

    Parameters
    ----------
    min: array_like
        Minutes to be converted.

    Returns
    -------
    seconds: float or array of floats
        Equivalent time in seconds.

    Notes
    -----
    Computes `` sec = min * 60 `` where 1 minute = 60 seconds .

    """
    return _npasanyarray(min)*60

def sec_to_h(sec):
    """
    Convert seconds to hours

    Parameters
    ----------
    seconds: array_like
        Seconds to be converted.

    Returns
    -------
    hours: float or array of floats
        Equivalent time in hours.

    Notes
    -----
    Computes `` h = sec / 3600 `` where 1h = 3600 sec .

    """
    return _npasanyarray(sec)/3600

def h_to_sec(h):
    """
    Convert hours to seconds

    Parameters
    ----------
    h: array_like
        Hours to be converted.

    Returns
    -------
    Seconds: float or array of floats
        Equivalent time in seconds.

    Notes
    -----
    Computes `` sec = 3600*h `` where 1h = 3600 sec.

    """
    return _npasanyarray(h)*3600

def sec_to_day(sec):
    """
    Convert seconds to day

    Parameters
    ----------
    sec: array_like
        Seconds to be converted.

    Returns
    -------
    day: float or array of floats
        Equivalent time in days.

    Notes
    -----
    Computes `` day = sec / 24*3600 `` where 1 day = 24*3600.

    """
    return _npasanyarray(sec)/24*3600

def day_to_sec(day):
    """
    Convert days to seconds

    Parameters
    ----------
    day: array_like
        Days to be converted.

    Returns
    -------
    seconds: float or array of floats
        Equivalent time in seconds.

    Notes
    -----
    Computes `` sec = day*24*3600 `` where 1 day = 24*3600.

    """
    return _npasanyarray(day)*24*3600

def sec_to_week(sec):
    """
    Convert seconds to weeks

    Parameters
    ----------
    sec: array_like
        seconds to be converted.

    Returns
    -------
    week: float or array of floats
        Equivalent time in weeks.

    Notes
    -----
    Computes `` week = sec/7*24*3600 `` where 1 week = 7*24*3600.

    """
    return _npasanyarray(sec)/7*24*3600

def week_to_sec(week):
    """
    Convert weeks to seconds

    Parameters
    ----------
    week: array_like
        Weeks to be converted.

    Returns
    -------
    Seconds: float or array of floats
        Equivalent time in weeks.

    Notes
    -----
    Computes `` sec = week*7*24*3600 `` where 1 week = 7*24*3600.

    """
    return _npasanyarray(week)*7*24*3600

def sec_to_year(year):
    """
    Convert seconds to years

    Parameters
    ----------
    sec: array_like
        Seconds to be converted.

    Returns
    -------
    year: float or array of floats
        Equivalent time in years.

    Notes
    -----
    Computes `` year =  sec / 365*24*3600 `` where 1 year = 365*24*3600 .

    """
    return _npasanyarray(year)/365*24*3600

def year_to_sec(year):
    """
    Convert years to seconds

    Parameters
    ----------
    year: array_like
        years to be converted.

    Returns
    -------
    seconds: float or array of floats
        Equivalent time in seconds.

    Notes
    -----
    Computes `` sec = year * 365*24*3600 `` where 1 year = 365*24*3600.

    """
    return _npasanyarray(sec)*365*24*3600

def min_to_h(min):
    """
    Convert minutes to hours

    Parameters
    ----------
    min: array_like
        minutes to be converted.

    Returns
    -------
    hours: float or array of floats
        Equivalent time in hours.

    Notes
    -----
    Computes `` h = min/60 `` where 1h = 60min .

    """
    return _npasanyarray(min)/60

def h_to_min(h):
    """
    Convert hours to minutes

    Parameters
    ----------
    h: array_like
        hours to be converted.

    Returns
    -------
    min: float or array of floats
        Equivalent time in minutes.

    Notes
    -----
    Computes `` min = h*60 `` where 1X = h = 60min.

    """
    return _npasanyarray(min)*60


def min_to_day(min):
    """
    Convert minutes to days

    Parameters
    ----------
    min: array_like
        minutes to be converted.

    Returns
    -------
    day: float or array of floats
        Equivalent time in days.

    Notes
    -----
    Computes `` day = min/24*60`` where 1 day = 24*60 min.

    """
    return _npasanyarray(min)/24*60

def day_to_min(day):
    """
    Convert days to minutes

    Parameters
    ----------
    day: array_like
        Days to be converted.

    Returns
    -------
    min: float or array of floats
        Equivalent time in minutes.

    Notes
    -----
    Computes `` min = day*24*60 `` where 1X = where 1 day = 24*60 min.

    """
    return _npasanyarray(day)*24*60

def min_to_week(min):
    """
    Convert minutes to weeks

    Parameters
    ----------
    min: array_like
        Minutes to be converted.

    Returns
    -------
    week: float or array of floats
        Equivalent time in weeks.

    Notes
    -----
    Computes `` week = min/7*24*60 `` where 1 week = 7*24*60 min.

    """
    return _npasanyarray(min)/7*24*60

def week_to_min(week):
    """
    Convert weeks to minutes

    Parameters
    ----------
    week: array_like
        Weeks to be converted.

    Returns
    -------
    min: float or array of floats
        Equivalent time in minutes.

    Notes
    -----
    Computes `` min = week*7*24*60 min `` where 1 week = 7*24*60 min .

    """
    return _npasanyarray(week)*7*24*60

def min_to_year(min):
    """
    Convert minutes to years

    Parameters
    ----------
    min: array_like
        Minutes to be converted.

    Returns
    -------
    year: float or array of floats
        Equivalent time in years.

    Notes
    -----
    Computes `` year = min/365*24*60 `` where 1 year = 365*24*60 .

    """
    return _npasanyarray(min)/365*24*60

def year_to_min(year):
    """
    Convert years to minutes

    Parameters
    ----------
    year: array_like
        Years to be converted.

    Returns
    -------
    min: float or array of floats
        Equivalent time in minutes.

    Notes
    -----
    Computes `` min = year*365*24*60 `` where 1 year = 365*24*60.

    """
    return _npasanyarray(year)*365*24*60

def h_to_day(h):
    """
    Convert hours to days.

    Parameters
    ----------
    h: array_like
        Hours to be converted.

    Returns
    -------
    days: float or array of floats
        Equivalent time in days.

    Notes
    -----
    Computes `` day = h/24`` where 1 day  = 24h .

    """
    return _npasanyarray(h)/24

def day_to_h(day):
    """
    Convert days to hours.

    Parameters
    ----------
    day: array_like
        Days to be converted.

    Returns
    -------
    h: float or array of floats
        Equivalent time in hours.

    Notes
    -----
    Computes `` h = 24*day `` where 1 day  = 24h.

    """
    return _npasanyarray(day)*24


def h_to_week(h):
    """
    Convert hours to weeks.

    Parameters
    ----------
    h: array_like
        Hours to be converted.

    Returns
    -------
    weeks: float or array of floats
        Equivalent time in weeks.

    Notes
    -----
    Computes `` week = h/24*7 `` where 1 week = 24*7h.

    """
    return _npasanyarray(h)/24*7

def week_to_h(week):
    """
    Convert weeks to hours

    Parameters
    ----------
    week: array_like
        Weeks to be converted.

    Returns
    -------
    h: float or array of floats
        Equivalent time in hours.

    Notes
    -----
    Computes `` h = week*24*7 `` where 1 week = 24*7h.

    """
    return _npasanyarray(week)*24*7

def h_to_year(h):
    """
    Convert hours to years

    Parameters
    ----------
    h: array_like
        Hours to be converted.

    Returns
    -------
    year: float or array of floats
        Equivalent time in years.

    Notes
    -----
    Computes `` year = h/365*24 `` where 1 year = 365*24h .

    """
    return _npasanyarray(h)/365*24

def year_to_h(h):
    """
    Convert years to hours

    Parameters
    ----------
    year: array_like
        Years to be converted.

    Returns
    -------
    h: float or array of floats
        Equivalent time in hours.

    Notes
    -----
    Computes `` h = year*365*24 `` where 1 year = 365*24h .

    """
    return _npasanyarray(h)*365*24

def day_to_week(day):
    """
    Convert days to weeks

    Parameters
    ----------
    day: array_like
        Days to be converted.

    Returns
    -------
    week: float or array of floats
        Equivalent time in weeks.

    Notes
    -----
    Computes `` week = day/7 `` where 1 week = 7 days .

    """
    return _npasanyarray(day)/7

def week_to_day(week):
    """
    Convert days to weeks

    Parameters
    ----------
    week: array_like
        Weeks to be converted.

    Returns
    -------
    day: float or array of floats
        Equivalent time in days.

    Notes
    -----
    Computes `` day = week*7 `` where 1 week = 7 days .

    """
    return _npasanyarray(day)*7