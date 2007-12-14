"Sparse Matrix Support"

from info import __doc__

from base import *
from compressed import *
from lil import *
from dok import *
from coo import *
from dia import *

from construct import *
from spfuncs import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import NumpyTest
test = NumpyTest().test
