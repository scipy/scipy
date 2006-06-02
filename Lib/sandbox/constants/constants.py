"""
Collection of physical constants and conversion factors.

Most constants are in SI units, so you can do
print '10 mile per minute is', 10*mile/minute, 'm/s or', 10*mile/(minute*knot), 'knots'

The list is not meant to be comprehensive, but just a convenient list for everyday use.
"""

"""
BasSw 2006
physical constants: imported from CODATA
unit conversion: see e.g. NIST special publication 811
Use at own risk: double-check values before calculating your Mars orbit-insertion burn.
"""

import math as _math
from codata import value as _cd

#mathematical constants
pi = _math.pi
golden = golden_ratio = (1 + _math.sqrt(5)) / 2

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

#physical constants

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
k = Bolzmann = _cd('Boltzmann constant')
sigma = Stefan_Bolzmann = _cd('Stefan-Boltzmann constant')
Wien = _cd('Wien displacement law constant')
Rydberg = _cd('Rydberg constant')

#weight in kg
lb = pound = 0.45359237 #definition
oz = ounce = pound / 16
stone = 14 * pound
long_ton = 2240 * pound
short_ton = 2000 * pound
metric_ton = 1000.0
m_e = electron_mass = _cd('electron mass')
m_p = proton_mass = _cd('proton mass')
m_n = neutron_mass = _cd('neutron mass')
m_u = u = atomic_mass = _cd('atomic mass constant')

#angle in rad
degree = pi / 180
arcmin = arcminute = degree / 60
arcsec = arcsecond = arcmin / 60

#time in second
minute = 60.0
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day
Julian_year = 365.25 * day

#length in meter
inch = 0.0254
foot = 12 * inch
yard = 3 * foot
mile = 1760 * yard
nm = nautical_mile = 1852.0
angstrom = 1e-10
micron = 1e-6
au = astronimical_unit = 149597870691.0
parsec = au / arcsec
lightyear = Julian_year * c

#pressure in Pascal
atm = atmosphere = _cd('standard atmosphere')
bar = 1e5 
torr = mmHg = atm / 760
psi = pound * g / (inch * inch)

#astronomy
#M_earth
#M_sun
#M_jup
#R_earth = earth_radius

#area in meter**2
hectare = 1e4
acre = 43560 * foot**2

#volume in meter**3
litre = liter = 1e-3
gallon = 231 * inch**3
pint = gallon / 8

#speed in meter per second
kmh = 1e3 / hour
mph = mile / hour
mach = speed_of_sound = 340.5 #approx value at 15 degrees in 1 atm. is this a common value?
knot = nm / hour

#energy in joule
eV = electron_volt = elementary_charge # * 1 Volt
calorie = 4.184 # or 4.1868?

#power in watt
hp = horsepower = 550 * foot * pound * g

#temperature conversions 
zero_Celsius = 273.15 #Kelvin

def C2K(C):
    """Convert Celcius to Kelvin"""
    return C + zero_Celsius

def K2C(K):
    """Convert Kelvin to Celcius"""
    return K - zero_Celsius

def F2C(F):
    """Convert Fahrenheit to Celcius"""
    return (F - 32) / 1.8

def C2F(C):
    """Convert Celcius to Fahrenheit"""
    return 1.8 * C + 32

def F2K(F):
    """Convert Fahrenheit to Kelvin"""
    return C2K(F2C(F))

def K2F(k):
    """Convert Kelvin to Fahrenheit"""
    return C2F(K2C(K))

#optics

def lambda2nu(lambda_):
    """Convert wavelength to optical frequency"""
    return c / lambda_

def nu2lambda(nu):
    """Convert optical frequency to wavelength"""
    return c / nu
