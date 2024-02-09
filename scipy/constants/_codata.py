"""
Fundamental Physical Constants
------------------------------

These constants are taken from CODATA Recommended Values of the Fundamental
Physical Constants 2018.

Object
------
physical_constants : dict
    A dictionary containing physical constants. Keys are the names of physical
    constants, values are tuples (value, units, precision).

Functions
---------
value(key):
    Returns the value of the physical constant(key).
unit(key):
    Returns the units of the physical constant(key).
precision(key):
    Returns the relative precision of the physical constant(key).
find(sub):
    Prints or returns list of keys containing the string sub, default is all.

Source
------
The values of the constants provided at this site are recommended for
international use by CODATA and are the latest available. Termed the "2018
CODATA recommended values," they are generally recognized worldwide for use in
all fields of science and technology. The values became available on 20 May
2019 and replaced the 2014 CODATA set. Also available is an introduction to the
constants for non-experts at

https://physics.nist.gov/cuu/Constants/introduction.html

References
----------
Theoretical and experimental publications relevant to the fundamental constants
and closely related precision measurements published since the mid 1980s, but
also including many older papers of particular interest, some of which date
back to the 1800s. To search the bibliography, visit

https://physics.nist.gov/cuu/Constants/

"""

# Compiled by Charles Harris, dated October 3, 2002
# updated to 2002 values by BasSw, 2006
# Updated to 2006 values by Vincent Davis June 2010
# Updated to 2014 values by Joseph Booker, 2015
# Updated to 2018 values by Jakob Jakobson, 2019

from __future__ import annotations

import warnings

from typing import Any

__all__ = ['physical_constants', 'value', 'unit', 'precision', 'find',
           'ConstantWarning']

"""
Source:  https://physics.nist.gov/cuu/Constants/

The values of the constants provided at this site are recommended for
international use by CODATA and are the latest available. Termed the "2018
CODATA recommended values," they are generally recognized worldwide for use in
all fields of science and technology. The values became available on 20 May
2019 and replaced the 2014 CODATA set.
"""

# https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2002.txt
with open('scipy/scipy/constants/codata_constants_2002.txt') as f:
    txt2002 = f.read()

# https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2006.txt
with open('scipy/scipy/constants/codata_constants_2006.txt') as f:
    txt2006 = f.read()

# https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2010.txt
with open('scipy/scipy/constants/codata_constants_2010.txt') as f:
    txt2010 = f.read()

# https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2014.txt
with open('scipy/scipy/constants/codata_constants_2014.txt') as f:
    txt2014 = f.read()

# https://physics.nist.gov/cuu/Constants/Table/allascii.txt
with open('scipy/scipy/constants/codata_constants_2018.txt') as f:
    txt2018 = f.read()

physical_constants: dict[str, tuple[float, str, float]] = {}


def parse_constants_2002to2014(d: str) -> dict[str, tuple[float, str, float]]:
    constants = {}
    for line in d.split('\n'):
        if line == "":
            continue
        if line[0] == " ":
            continue
        if line[0] == "-":
            continue        
        name = line[:55].rstrip()
        val = float(line[55:77].replace(' ', '').replace('...', ''))
        uncert = float(line[77:99].replace(' ', '').replace('(exact)', '0'))
        units = line[99:].rstrip()
        constants[name] = (val, units, uncert)
    return constants


def parse_constants_2018toXXXX(d: str) -> dict[str, tuple[float, str, float]]:
    constants = {}
    for line in d.split('\n'):
        if line == "":
            continue
        if line[0] == " ":
            continue
        if line[0] == "-":
            continue        
        name = line[:60].rstrip()
        val = float(line[60:85].replace(' ', '').replace('...', ''))
        uncert = float(line[85:110].replace(' ', '').replace('(exact)', '0'))
        units = line[110:].rstrip()
        constants[name] = (val, units, uncert)
    return constants

_physical_constants_2002 = parse_constants_2002to2014(txt2002)
_physical_constants_2006 = parse_constants_2002to2014(txt2006)
_physical_constants_2010 = parse_constants_2002to2014(txt2010)
_physical_constants_2014 = parse_constants_2002to2014(txt2014)
_physical_constants_2018 = parse_constants_2018toXXXX(txt2018)

physical_constants.update(_physical_constants_2002)
physical_constants.update(_physical_constants_2006)
physical_constants.update(_physical_constants_2010)
physical_constants.update(_physical_constants_2014)
physical_constants.update(_physical_constants_2018)
_current_constants = _physical_constants_2018
_current_codata = "CODATA 2018"

# check obsolete values
_obsolete_constants = {}
for k in physical_constants:
    if k not in _current_constants:
        _obsolete_constants[k] = True

# generate some additional aliases
_aliases = {}
for k in _physical_constants_2002:
    if 'magn.' in k:
        _aliases[k] = k.replace('magn.', 'mag.')
for k in _physical_constants_2006:
    if 'momentum' in k:
        _aliases[k] = k.replace('momentum', 'mom.um')
for k in _physical_constants_2018:
    if 'momentum' in k:
        _aliases[k] = k.replace('momentum', 'mom.um')

# CODATA 2018: renamed and no longer exact; use as aliases
_aliases['mag. constant'] = 'vacuum mag. permeability'
_aliases['electric constant'] = 'vacuum electric permittivity'


class ConstantWarning(DeprecationWarning):
    """Accessing a constant no longer in current CODATA data set"""
    pass


def _check_obsolete(key: str) -> None:
    if key in _obsolete_constants and key not in _aliases:
        warnings.warn(f"Constant '{key}' is not in current {_current_codata} data set",
                      ConstantWarning, stacklevel=3)


def value(key: str) -> float:
    """
    Value in physical_constants indexed by key

    Parameters
    ----------
    key : Python string
        Key in dictionary `physical_constants`

    Returns
    -------
    value : float
        Value in `physical_constants` corresponding to `key`

    Examples
    --------
    >>> from scipy import constants
    >>> constants.value('elementary charge')
    1.602176634e-19

    """
    _check_obsolete(key)
    return physical_constants[key][0]


def unit(key: str) -> str:
    """
    Unit in physical_constants indexed by key

    Parameters
    ----------
    key : Python string
        Key in dictionary `physical_constants`

    Returns
    -------
    unit : Python string
        Unit in `physical_constants` corresponding to `key`

    Examples
    --------
    >>> from scipy import constants
    >>> constants.unit('proton mass')
    'kg'

    """
    _check_obsolete(key)
    return physical_constants[key][1]


def precision(key: str) -> float:
    """
    Relative precision in physical_constants indexed by key

    Parameters
    ----------
    key : Python string
        Key in dictionary `physical_constants`

    Returns
    -------
    prec : float
        Relative precision in `physical_constants` corresponding to `key`

    Examples
    --------
    >>> from scipy import constants
    >>> constants.precision('proton mass')
    5.1e-37

    """
    _check_obsolete(key)
    return physical_constants[key][2] / physical_constants[key][0]


def find(sub: str | None = None, disp: bool = False) -> Any:
    """
    Return list of physical_constant keys containing a given string.

    Parameters
    ----------
    sub : str
        Sub-string to search keys for. By default, return all keys.
    disp : bool
        If True, print the keys that are found and return None.
        Otherwise, return the list of keys without printing anything.

    Returns
    -------
    keys : list or None
        If `disp` is False, the list of keys is returned.
        Otherwise, None is returned.

    Examples
    --------
    >>> from scipy.constants import find, physical_constants

    Which keys in the ``physical_constants`` dictionary contain 'boltzmann'?

    >>> find('boltzmann')
    ['Boltzmann constant',
     'Boltzmann constant in Hz/K',
     'Boltzmann constant in eV/K',
     'Boltzmann constant in inverse meter per kelvin',
     'Stefan-Boltzmann constant']

    Get the constant called 'Boltzmann constant in Hz/K':

    >>> physical_constants['Boltzmann constant in Hz/K']
    (20836619120.0, 'Hz K^-1', 0.0)

    Find constants with 'radius' in the key:

    >>> find('radius')
    ['Bohr radius',
     'classical electron radius',
     'deuteron rms charge radius',
     'proton rms charge radius']
    >>> physical_constants['classical electron radius']
    (2.8179403262e-15, 'm', 1.3e-24)

    """
    if sub is None:
        result = list(_current_constants.keys())
    else:
        result = [key for key in _current_constants
                  if sub.lower() in key.lower()]

    result.sort()
    if disp:
        for key in result:
            print(key)
        return
    else:
        return result


c = value('speed of light in vacuum')
mu0 = value('vacuum mag. permeability')
epsilon0 = value('vacuum electric permittivity')

# Table is lacking some digits for exact values: calculate from definition
exact_values = {
    'joule-kilogram relationship': (1 / (c * c), 'kg', 0.0),
    'kilogram-joule relationship': (c * c, 'J', 0.0),
    'hertz-inverse meter relationship': (1 / c, 'm^-1', 0.0),
}

# sanity check
for key in exact_values:
    val = physical_constants[key][0]
    if abs(exact_values[key][0] - val) / val > 1e-9:
        raise ValueError("Constants.codata: exact values too far off.")
    if exact_values[key][2] == 0 and physical_constants[key][2] != 0:
        raise ValueError("Constants.codata: value not exact")

physical_constants.update(exact_values)

_tested_keys = ['natural unit of velocity',
                'natural unit of action',
                'natural unit of action in eV s',
                'natural unit of mass',
                'natural unit of energy',
                'natural unit of energy in MeV',
                'natural unit of mom.um',
                'natural unit of mom.um in MeV/c',
                'natural unit of length',
                'natural unit of time']

# finally, insert aliases for values
for k, v in list(_aliases.items()):
    if v in _current_constants or v in _tested_keys:
        physical_constants[k] = physical_constants[v]
    else:
        del _aliases[k]
