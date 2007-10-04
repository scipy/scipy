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

import unittest
import md5
import tempfile
import shutil
from os.path import join

import numpy as np
from scipy.io.nifti import NiftiImage


def md5sum(filename):
    """ Generate MD5 hash string.
    """
    file = open( filename )
    sum = md5.new()
    while True:
        data = file.read()
        if not data:
            break
        sum.update(data)
    return sum.hexdigest()


class TestNiftiImage(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp('pynifti_test')

    def tearDown(self):
        shutil.rmtree(self.workdir)

    def test_idempotent_load_save_cycle(self):
        """ check if file is unchanged by load/save cycle.
        """
        md5_orig = md5sum('data/example4d.nii.gz')
        nimg = NiftiImage('data/example4d.nii.gz')
        nimg.save( join( self.workdir, 'iotest.nii.gz') )
        md5_io =  md5sum( join( self.workdir, 'iotest.nii.gz') )

        self.failUnlessEqual(md5_orig, md5_io)

    def test_qform_setting(self):
        nimg = NiftiImage('data/example4d.nii.gz')
        # 4x4 identity matrix
        ident = np.identity(4)
        self.failIf( (nimg.qform == ident).all() )

        # assign new qform
        nimg.qform = ident
        self.failUnless( (nimg.qform == ident).all() )

        # test save/load cycle
        nimg.save( join( self.workdir, 'qformtest.nii.gz') )
        nimg2 = NiftiImage( join( self.workdir, 'qformtest.nii.gz') )

        self.failUnless( (nimg.qform == nimg2.qform).all() )


if __name__ == '__main__':
    NumpyTest.run()

