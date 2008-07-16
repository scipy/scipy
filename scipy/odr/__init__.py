#
# odr - Orthogonal Distance Regression
#

from info import __doc__

__version__ = '0.7'
__author__ = 'Robert Kern <robert.kern@gmail.com>'
__date__ = '2006-09-21'

import odrpack
from odrpack import odr         ,\
                    odr_error   ,\
                    odr_stop    ,\
                    Data        ,\
                    RealData    ,\
                    Model       ,\
                    Output      ,\
                    ODR

__all__ = ['odr', 'odr_error', 'odr_stop', 'Data', 'RealData', 'Model',
           'Output', 'ODR', 'odrpack']

from numpy.testing import Tester
test = Tester().test
#### EOF #######################################################################
