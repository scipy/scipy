from .codata import value
import numpy as np

# auxiliary functions
def convert_temperature(val, old_scale, new_scale):
    """
    Convert from a temperature scale to another one among Celsius, Kelvin,
    Fahrenheit, and Rankine scales.

    Parameters
    ----------
    val : array_like
        Value(s) of the temperature(s) to be converted expressed in the
        original scale.

    old_scale: str
        Specifies as a string the original scale from which the temperature
        value(s) will be converted. Supported scales are Celsius ('Celsius',
        'celsius', 'C' or 'c'), Kelvin ('Kelvin', 'kelvin', 'K', 'k'),
        Fahrenheit ('Fahrenheit', 'fahrenheit', 'F' or 'f'), and Rankine
        ('Rankine', 'rankine', 'R', 'r').

    new_scale: str
        Specifies as a string the new scale to which the temperature
        value(s) will be converted. Supported scales are Celsius ('Celsius',
        'celsius', 'C' or 'c'), Kelvin ('Kelvin', 'kelvin', 'K', 'k'),
        Fahrenheit ('Fahrenheit', 'fahrenheit', 'F' or 'f'), and Rankine
        ('Rankine', 'rankine', 'R', 'r').

    Returns
    -------
    res : float or array of floats
        Value(s) of the converted temperature(s) expressed in the new scale.

    Notes
    -----
    .. versionadded:: 0.18.0

    Examples
    --------
    >>> from scipy.constants import convert_temperature
    >>> convert_temperature(np.array([-40, 40.0]), 'Celsius', 'Kelvin')
    array([ 233.15,  313.15])

    """
    # Convert from `old_scale` to Kelvin
    if old_scale.lower() in ['celsius', 'c']:
        tempo = _np.asanyarray(val) + 273.15
    elif old_scale.lower() in ['kelvin', 'k']:
        tempo = _np.asanyarray(val)
    elif old_scale.lower() in ['fahrenheit', 'f']:
        tempo = (_np.asanyarray(val) - 32) * 5 / 9 + 273.15
    elif old_scale.lower() in ['rankine', 'r']:
        tempo = _np.asanyarray(val) * 5 / 9
    else:
        raise NotImplementedError("%s scale is unsupported: supported scales "
                                  "are Celsius, Kelvin, Fahrenheit, and "
                                  "Rankine" % old_scale)
    # and from Kelvin to `new_scale`.
    if new_scale.lower() in ['celsius', 'c']:
        res = tempo - 273.15
    elif new_scale.lower() in ['kelvin', 'k']:
        res = tempo
    elif new_scale.lower() in ['fahrenheit', 'f']:
        res = (tempo - 273.15) * 9 / 5 + 32
    elif new_scale.lower() in ['rankine', 'r']:
        res = tempo * 9 / 5
    else:
        raise NotImplementedError("'%s' scale is unsupported: supported "
                                  "scales are 'Celsius', 'Kelvin', "
                                  "'Fahrenheit', and 'Rankine'" % new_scale)

    return res


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
    >>> import scipy.constants as sc
    >>> sc.lambda2nu(np.array((1, sc.value('speed of light in vacuum'))))
    array([  2.99792458e+08,   1.00000000e+00])

    """
    return np.asanyarray(value('speed of light in vacuum')) / lambda_


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
    >>> import scipy.constants as sc
    >>> sc.nu2lambda(np.array((1, sc.value('speed of light in vacuum'))))
    array([  2.99792458e+08,   1.00000000e+00])

    """
    return value('speed of light in vacuum') / np.asanyarray(nu)
