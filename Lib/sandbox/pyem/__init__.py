#! /usr/bin/env python
# Last Change: Sun Jul 22 11:00 AM 2007 J

raise ImportError(
"""pyem has been moved to scikits and renamed to em. Please install
scikits.learn instead, and change your import to the following:

from scickits.learn.machine import em.""")

from info import __doc__

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
#from online_em import OnGMM as _OnGMM
#import examples as _examples

__all__ = filter(lambda s:not s.startswith('_'), dir())

from numpy.testing import NumpyTest
test = NumpyTest().test

