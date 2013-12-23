"""Utility to compare (Numpy) version strings.

The NumpyVersion class allows properly comparing numpy version strings.
The LooseVersion and StrictVersion classes that distutils provides don't
work; they don't recognize anything like alpha/beta/rc/dev versions.

"""

import re

from scipy.lib.six import string_types


__all__ = ['NumpyVersion']


class NumpyVersion():
    """Parse and compare numpy version strings.

    Numpy has the following versioning scheme (numbers given are examples; they
    can be >9) in principle):

    - Released version: '1.8.0', '1.8.1', etc.
    - Alpha: '1.8.0a1', '1.8.0a2', etc.
    - Beta: '1.8.0b1', '1.8.0b2', etc.
    - Release candidates: '1.8.0rc1', '1.8.0rc2', etc.
    - Development versions: '1.8.0.dev-f1234afa' (git commit hash appended)
    - Development versions after a1: '1.8.0a1.dev-f1234afa',
                                     '1.8.0b2.dev-f1234afa',
                                     '1.8.1rc1.dev-f1234afa', etc.
    - Development versions (no git hash available): '1.8.0.dev-Unknown'

    Comparing needs to be done against a valid version string or other
    `NumpyVersion` instance.

    Parameters
    ----------
    vstring : str
        Numpy version string (``np.__version__``).

    Notes
    -----
    All dev versions of the same (pre-)release compare equal.

    Examples
    --------
    >>> from scipy.lib._version import NumpyVersion
    >>> if NumpyVersion(np.__version__) < '1.7.0'):
    ...     print('skip')
    skip

    >>> NumpyVersion('1.7')  # raises ValueError, add ".0"

    """
    def __init__(self, vstring):
        self.vstring = vstring
        if not re.match(r'\d[.]\d[.]\d', vstring):
            raise ValueError("Not a valid numpy version string")

        self.version = vstring[:5]
        if len(vstring) == 5:
            self.pre_release = 'z_final'  # compare > (a, b, rc)
        else:
            alpha = re.match(r'a\d', vstring[5:])
            beta = re.match(r'b\d', vstring[5:])
            rc = re.match(r'rc\d', vstring[5:])
            pre_rel = [m for m in [alpha, beta, rc] if m is not None]
            if pre_rel:
                self.pre_release = pre_rel[0].group()
            else:
                self.pre_release = None

        self.is_devversion = bool(re.match(r'.dev-', vstring))


    def __cmp__(self, other):
        if not isinstance(other, (string_types, NumpyVersion)):
            raise ValueError("Invalid object to compare with NumpyVersion.")

        if isinstance(other, string_types):
            other = NumpyVersion(other)

        vercmp = cmp(self.version, other.version)
        if not vercmp == 0:
            return vercmp

        # Same x.y.z version, check for alpha/beta/rc
        vercmp = cmp(self.pre_release, other.pre_release)
        if not vercmp == 0:
            return vercmp

        # Same version and same pre-release, check if dev version
        if self.is_devversion and not other.is_devversion:
            return -1
        elif self.is_devversion and other.is_devversion:
            return 0
        elif not self.is_devversion and other.is_devversion:
            return 1
        else:
            return 0
