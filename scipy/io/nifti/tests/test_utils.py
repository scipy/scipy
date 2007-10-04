### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#
#    Unit tests for PyNIfTI file io
#
#    Copyright (C) 2007 by
#    Michael Hanke <michael.hanke@gmail.com>
#
#    This is free software; you can redistribute it and/or
#    modify it under the terms of the MIT License.
#
#    This package is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the COPYING
#    file that comes with this package for more details.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

import scipy.io.nifti.utils
import numpy
import unittest


class UtilsTests(unittest.TestCase):
    def testZScoring(self):
        # dataset: mean=2, std=1
        data = numpy.array( (0,1,3,4,2,2,3,1,1,3,3,1,2,2,2,2) )
        self.failUnlessEqual( data.mean(), 2.0 )
        self.failUnlessEqual( data.std(), 1.0 )


def suite():
    return unittest.makeSuite(UtilsTests)


if __name__ == '__main__':
    unittest.main()

