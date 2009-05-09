# lineshape functors
# -*- coding: iso-8859-1 -*-
#
# Copyright (C) 2002 Jochen Küpper <jochen@jochen-kuepper.de>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


__doc__     = """Lineshape functors.

The objects defined in this module can be used to calculate numarrays containing
common lineshape-profiles.

For the *Profile classes the profile is only evaluated once for each object and
then reused for each call. If you only use a profile once and you are concerned
about memory consumption you could call the underlying functions directly.
"""

__author__  = "Jochen Küpper <jochen@jochen-kuepper.de>"
__date__    = "$Date: 2007/03/14 16:35:57 $"[7:-11]
__version__ = "$Revision: 1.1 $"[11:-2]

from convolve._lineshape import *

class Profile(object):
    """An base object to provide a convolution kernel."""

    def __init__(self, x, w, x0=0.0):
        # call init for all superclasses
        super(Profile, self).__init__(x, w, x0)
        self._recalculate(x, w, x0)

    def __call__(self):
        return self._kernel

    def _recalculate(self, x, w, x0):
        self._kernel = self._profile(x, w, x0)


class GaussProfile(Profile):
    """An object for Gauss-folding."""

    def __init__(self, x, w, x0=0.0):
        self._profile = gauss
        # call init for all superclasses
        super(GaussProfile, self).__init__(x, w, x0)


class LorentzProfile(Profile):
    """An object for Lorentz-folding."""

    def __init__(self, x, w, x0=0.0):
        self._profile = lorentz
        # call init for all superclasses
        super(LorentzProfile, self).__init__(x, w, x0)



class VoigtProfile(Profile):
    """An object for Voigt-folding.

    The constructor takes the following parameter:
      |x|  Scalar or numarray with values to calculate profile at.
      |w|  Tuple of Gaussian and Lorentzian linewidth contribution
      |x0| Center frequency
    """

    def __init__(self, x, w, x0=0.0):
        self._profile = voigt
        # call init for all superclasses
        super(VoigtProfile, self).__init__(x, w, x0)



## Local Variables:
## mode: python
## mode: auto-fill
## fill-column: 80
## End:
