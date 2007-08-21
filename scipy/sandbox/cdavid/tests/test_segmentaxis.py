#! /usr/bin/env python
# Last Change: Fri Dec 15 10:00 PM 2006 J

from numpy.testing import *

import numpy as N

set_package_path()
from cdavid.segmentaxis import segment_axis
restore_path()

# #Optional:
# set_local_path()
# # import modules that are located in the same directory as this file.
# restore_path()

class test_segment(NumpyTestCase):
    def check_simple(self):
        assert_equal(segment_axis(N.arange(6),length=3,overlap=0),
                         N.array([[0,1,2],[3,4,5]]))

        assert_equal(segment_axis(N.arange(7),length=3,overlap=1),
                         N.array([[0,1,2],[2,3,4],[4,5,6]]))

        assert_equal(segment_axis(N.arange(7),length=3,overlap=2),
                         N.array([[0,1,2],[1,2,3],[2,3,4],[3,4,5],[4,5,6]]))

    def check_error_checking(self):
        self.assertRaises(ValueError,
                lambda: segment_axis(N.arange(7),length=3,overlap=-1))
        self.assertRaises(ValueError,
                lambda: segment_axis(N.arange(7),length=0,overlap=0))
        self.assertRaises(ValueError,
                lambda: segment_axis(N.arange(7),length=3,overlap=3))
        self.assertRaises(ValueError,
                lambda: segment_axis(N.arange(7),length=8,overlap=3))

    def check_ending(self):
        assert_equal(segment_axis(N.arange(6),length=3,overlap=1,end='cut'),
                         N.array([[0,1,2],[2,3,4]]))
        assert_equal(segment_axis(N.arange(6),length=3,overlap=1,end='wrap'),
                         N.array([[0,1,2],[2,3,4],[4,5,0]]))
        assert_equal(segment_axis(N.arange(6),length=3,overlap=1,end='pad',endvalue=-17),
                         N.array([[0,1,2],[2,3,4],[4,5,-17]]))

    def check_multidimensional(self):
        
        assert_equal(segment_axis(N.ones((2,3,4,5,6)),axis=3,length=3,overlap=1).shape,
                     (2,3,4,2,3,6))

        assert_equal(segment_axis(N.ones((2,5,4,3,6)).swapaxes(1,3),axis=3,length=3,overlap=1).shape,
                     (2,3,4,2,3,6))

        assert_equal(segment_axis(N.ones((2,3,4,5,6)),axis=2,length=3,overlap=1,end='cut').shape,
                     (2,3,1,3,5,6))

        assert_equal(segment_axis(N.ones((2,3,4,5,6)),axis=2,length=3,overlap=1,end='wrap').shape,
                     (2,3,2,3,5,6))

        assert_equal(segment_axis(N.ones((2,3,4,5,6)),axis=2,length=3,overlap=1,end='pad').shape,
                     (2,3,2,3,5,6))

if __name__=='__main__':
    NumpyTest().run()
