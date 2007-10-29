"""
Wrapper to lowess and stl routines.

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date$
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'


import os

import numpy
from numpy import fromiter
from numpy import bool_, float_

#import maskedarray as MA

import pyloess
from pyloess import loess

com_example = [
"""
# Get some example data ...................................
dfile = open(os.path.join('tests','madeup_data'), 'r')
dfile.readline()
x = fromiter((float(v) for v in dfile.readline().rstrip().split()),
             float_).reshape(-1,2)
dfile.readline()
y = fromiter((float(v) for v in dfile.readline().rstrip().split()),
             float_)
""",

"""
# Get some additional info for prediction .................
newdata1 = numpy.array([[-2.5, 0.0, 2.5], [0., 0., 0.]])
newdata2 = numpy.array([[-0.5, 0.5], [0., 0.]])
""",

"""
# Create a new loess object ...............................
madeup = loess(x,y)
# ... and prints the parameters
print madeup.model,'\\n', madeup.control
""",

"""
# Modify some of the model parameters .....................
madeup.model.update(span=0.8, normalize=False)
print madeup.model
"""
]



if 1:
    for com in com_example:
        print com
        exec(com)
