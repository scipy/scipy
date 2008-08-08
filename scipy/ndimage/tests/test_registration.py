import math
import numpy as np
import scipy.ndimage._registration as reg
from numpy.testing import *

def load_desc():
    # this is for a 256x256x90 volume with 0.9375 x 0.9375 * 1.5 mm voxel sizes 
    rows   = 256
    cols   = 256
    layers = 90
    xsamp  = 0.9375
    ysamp  = 0.9375
    zsamp  = 1.5
    desc = {'rows' : rows, 'cols' : cols, 'layers' : layers, 
            'sample_x' : xsamp, 'sample_y' : ysamp, 'sample_z' : zsamp}
    return desc

def build_volume(imagedesc, S=[1500.0, 2500.0, 1000.0]):

    """
    build a 3D Gaussian volume. user passes in image dims in imagedesc
    the sigma for each axis is S[3] where 0=z, 1=y, 2=x

    volume3D = build_test_volume(imagedesc, S)

    Parameters 
    ----------
    imagedesc : {dictionary}
        volume dimensions and sampling

    S : {tuple}
        the Gaussian sigma for Z, Y and X

    Returns 
    -------

    volume3D : {nd_array}
        the 3D volume for testing

    """
    layers = imagedesc['layers']
    rows   = imagedesc['rows']
    cols   = imagedesc['cols']

    L = layers/2
    R = rows/2
    C = cols/2

    # build coordinates for 3D Gaussian volume
    # coordinates are centered at (0, 0, 0)
    [a, b, c] = np.mgrid[-L:L, -R:R, -C:C]

    sigma    = np.array([S[0], S[1], S[2]])
    aa       = (np.square(a))/sigma[0]
    bb       = (np.square(b))/sigma[1]
    cc       = (np.square(c))/sigma[2]
    volume3D = (255.0*np.exp(-(aa + bb + cc))).astype(np.uint8)

    return volume3D

 # self.failUnless(diff(output, tcov) < eps)
class TestRegistration(TestCase):

    def test_affine_matrix_build_1(self):
        "test_affine_matrix_build_1"
        P = np.zeros(6)
        M = reg.build_rotate_matrix(P)
        E = np.eye(4)
        match = (E==M).all()
        assert_equal(match, True)
        return

    def test_affine_matrix_build_2(self):
        "test_affine_matrix_build_2"
        P = np.zeros(6)
        P[0] = 1.0
        M = reg.build_rotate_matrix(P)
        E = np.array([
                     [ 1. ,  0.        , 0.        , 0. ],
                     [ 0. ,  0.9998477 , 0.01745241, 0. ],
                     [ 0. , -0.01745241, 0.9998477 , 0. ],
                     [ 0. ,  0.        , 0.        , 1. ]
                     ])
        assert_array_almost_equal(E, M, decimal=6)
        return

    def test_affine_matrix_build_3(self):
        "test_affine_matrix_build_3"
        P = np.zeros(6)
        P[0] = 1.0
        P[1] = 1.0
        P[2] = 1.0
        M = reg.build_rotate_matrix(P)
        E = np.array([
                     [ 0.99969541,  0.01744975,  0.01745241,  0. ],
                     [-0.01775429,  0.9996901 ,  0.01744975,  0. ],
                     [-0.0171425 , -0.01775429,  0.99969541,  0. ],
                     [ 0.        ,  0.        ,  0.        ,  1. ]
                     ])
        assert_array_almost_equal(E, M, decimal=6)
        return

    @dec.slow
    def test_autoalign_histogram_1(self):
        "test_autoalign_histogram_1"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost, joint = reg.check_alignment(gvol, mat, gvol, mat, ret_histo=1, lite=1)
        # confirm that with lite=1 only have non-zero on the main diagonal
        j_diag = joint.diagonal()
        Z = np.diag(j_diag)
        W = joint - Z
        # clip the near-zero fuzz
        W[abs(W) < 1e-10] = 0.0
        assert_equal(W.max(), 0.0)
        return

    @dec.slow
    def test_autoalign_histogram_2(self):
        "test_autoalign_histogram_2"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost, joint = reg.check_alignment(gvol, mat, gvol, mat, ret_histo=1, lite=0)
        # confirm that with lite=0 DO have non-zero on the main diagonal
        j_diag = joint.diagonal()
        Z = np.diag(j_diag)
        W = joint - Z
        # clip the near-zero fuzz
        W[abs(W) < 1e-10] = 0.0
        s = (W.max() > 0.0)
        # make sure there are off-diagonal values
        assert_equal(s, True)
        return

    @dec.slow
    def test_autoalign_ncc_value_1(self):
        "test_autoalign_ncc_value_1"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost = reg.check_alignment(gvol, mat, gvol, mat, method='ncc', lite=1)
        # confirm the cross correlation is near 1.0 
        t = abs(cost) + 0.0001
        s = (t >= 1.0)
        assert_equal(s, True)
        return

    @dec.slow
    def test_autoalign_ncc_value_2(self):
        "test_autoalign_ncc_value_2"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost = reg.check_alignment(gvol, mat, gvol, mat, method='ncc', lite=0)
        # confirm the cross correlation is near 1.0 
        t = abs(cost) + 0.0001
        s = (t >= 1.0)
        assert_equal(s, True)
        return

    @dec.slow
    def test_autoalign_nmi_value_1(self):
        "test_autoalign_nmi_value_1"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost = reg.check_alignment(gvol, mat, gvol, mat, method='nmi', lite=1)
        # confirm the normalized mutual information is near -2.0 
        assert_almost_equal(cost, -2.0, decimal=6)
        return

    @dec.slow
    def test_autoalign_nmi_value_2(self):
        "test_autoalign_nmi_value_2"
        desc = load_desc()
        gvol = build_volume(desc)
        mat  = np.eye(4)
        cost = reg.check_alignment(gvol, mat, gvol, mat, method='nmi', lite=0)
        # confirm the normalized mutual information is near -2.0 
        assert_almost_equal(cost, -1.7973048186515352, decimal=6)
        return



if __name__ == "__main__":
    #nose.runmodule()
    run_module_suite()



