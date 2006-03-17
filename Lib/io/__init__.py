
#
# io - Data input and output
#

from info import __doc__

from numpyio import packbits, unpackbits, bswap, fread, fwrite, \
     convert_objectarray
from mio import *
from array_import import *
from data_store import *
from pickler import *

from mmio import mminfo,mmread,mmwrite

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import ScipyTest
test = ScipyTest().test
