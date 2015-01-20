''' Byteorder utilities for system - numpy byteorder encoding

Converts a variety of string codes for little endian, big endian,
native byte order and swapped byte order to explicit numpy endian
codes - one of '<' (little endian) or '>' (big endian)

'''
from __future__ import division, print_function, absolute_import

import sys

sys_is_le = sys.byteorder == "little"
native_code = "<" if sys_is_le else ">"
swapped_code = ">" if sys_is_le else "<"

aliases = {"little": "<", "<": "<", "l": "<", "le": "<",
           "big": ">", ">": ">", "b": ">", "be": ">",
           None: native_code, "native": native_code, "=": native_code,
           "swapped": swapped_code, "s": swapped_code}


def to_numpy_code(code):
    """
    Convert various order codings to numpy format.

    Parameters
    ----------
    code : str
        The code to convert. It is converted to lower case before parsing.
        Legal values are:
        'little', 'big', 'l', 'b', 'le', 'be', '<', '>', 'native', '=',
        'swapped', 's'.

    Returns
    -------
    out_code : {'<', '>'}
        Here '<' is the numpy dtype code for little endian,
        and '>' is the code for big endian.

    Examples
    --------
    >>> import sys
    >>> sys_is_le == (sys.byteorder == 'little')
    True
    >>> to_numpy_code('big')
    '>'
    >>> to_numpy_code('little')
    '<'
    >>> nc = to_numpy_code('native')
    >>> nc == '<' if sys_is_le else nc == '>'
    True
    >>> sc = to_numpy_code('swapped')
    >>> sc == '>' if sys_is_le else sc == '<'
    True

    """
    try:
        return aliases[code.lower() if code else code]
    except KeyError:
        raise ValueError(
            'We cannot handle byte order %s' % code)
