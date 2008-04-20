import math
import numpy as NP
import scipy.ndimage._segmenter as seg
from scipy.testing import *

def run_sobel():
    img    = seg.build_test_discs()
    filter = seg.build_2d_kernel(hiFilterCutoff=60.0)
    fslice = seg.pre_filter(img, filter, low_threshold=0, high_threshold=255)
    sobel_edge_image, sobel_stats = seg.sobel_image(fslice)
    sobel_edge = seg.sobel_edges(sobel_edge_image, sobel_stats, sobel_threshold=0.8)
    label_sobel, sobel_groups = seg.get_blobs(sobel_edge)
    ROI      = seg.get_blob_regions(label_sobel, sobel_groups)
    measures = seg.get_all_bounding_boxes(ROI)
    thin_kernel = seg.build_morpho_thin_masks()
    sobel_edges = seg.mat_filter(label_sobel, thin_kernel, ROI)
    seg.get_voxel_measures(label_sobel, img, ROI)
    means = ROI[:]['voxelMean']
    return measures, means

def run_canny():
    img          = seg.build_test_discs()
    filter       = seg.build_2d_kernel(hiFilterCutoff=60.0)
    fslice       = seg.pre_filter(img, filter, low_threshold=0, high_threshold=255)
    canny_kernel = seg.build_d_gauss_kernel()
    horz, vert, imean = seg.canny_filter(fslice, canny_kernel)
    mag, canny_stats  = seg.canny_nonmax_supress(horz, vert, imean)
    canny_edge        = seg.canny_hysteresis(mag, canny_stats)
    label_canny, canny_groups = seg.get_blobs(canny_edge)
    ROI      = seg.get_blob_regions(label_canny, canny_groups)
    measures = seg.get_all_bounding_boxes(ROI)
    seg.get_voxel_measures(label_canny, img, ROI)
    means = ROI[:]['voxelMean']
    return measures, means

def run_texture1():
    filter      = seg.build_2d_kernel(hiFilterCutoff=60.0)
    img         = seg.build_test_discs()
    disc_mask   = seg.pre_filter(img, filter, low_threshold=50, high_threshold=255,
                                 conv_binary=1)
    label_disc_mask, disc_mask_groups = seg.get_blobs(disc_mask)
    disc_ROI    = seg.get_blob_regions(label_disc_mask, disc_mask_groups)
    laws_kernel = seg.build_laws_kernel() 
    impulse     = seg.build_test_impulses()
    calib       = seg.laws_texture_filter(impulse, label_disc_mask, laws_kernel,
                                          ROI=disc_ROI, verbose=1)
    kernels = calib[0]
    x = laws_kernel['coefficients'][0]
    m = NP.outer(x, x)
    m = m * 2
    laws_LL = kernels[0, 50-4:50+3, 50-3:50+4]
    return m, laws_LL 


def run_texture2():
    filter = seg.build_2d_kernel(hiFilterCutoff=60.0)
    img = seg.build_test_unit_discs()
    disc = seg.pre_filter(img, filter, low_threshold=50, high_threshold=255)
    disc_mask = seg.pre_filter(img, filter, low_threshold=50, high_threshold=255,
                               conv_binary=1)
    label_disc_mask, disc_mask_groups = seg.get_blobs(disc_mask)
    disc_ROI = seg.get_blob_regions(label_disc_mask, disc_mask_groups)
    laws_kernel = seg.build_laws_kernel() 
    texture_img = seg.build_test_texture_discs()
    seg.laws_texture_filter(texture_img, label_disc_mask, laws_kernel, ROI=disc_ROI,
                            mean_feature=1, verbose=0)
    tem = disc_ROI['TEM']
    return tem 

class TestSegment(TestCase):
    def test_sobel(self):
        # generate 4 discs, find the bounding boxes and 
        # confirm the bounding boxes are at the true position
        measures, voxel_means = run_sobel()
        number = measures.size
        _shortstruct = NP.dtype([('left', 'i'),
                                 ('right', 'i'),
                                 ('top', 'i'),
                                 ('bottom', 'i')])

        assert_equal(number, 4)
        # load the ground truth
        truth = NP.zeros(number, dtype=_shortstruct)
        truth[0] = (76, 179, 179, 77) 
        truth[1] = (332, 435, 179, 77)
        truth[2] = (76, 179, 435, 333)
        truth[3] = (332, 435, 435, 333)
        match = (truth==measures).all()
        assert_equal(match, True)
        # load the ground truth for the bounding box test image mean value
        voxel_truth = NP.zeros(number, dtype=NP.float64)
        voxel_truth = (80.0, 90.0, 100.0, 110.0)
        match = (voxel_truth==voxel_means).all()
        assert_equal(match, True)

        return

    def test_canny(self):
        # generate 4 discs, find the bounding boxes and 
        # confirm the bounding boxes are at the true position
        measures, voxel_means = run_canny()
        number = measures.size
        _shortstruct = NP.dtype([('left', 'i'),
                                 ('right', 'i'),
                                 ('top', 'i'),
                                 ('bottom', 'i')])

        assert_equal(number, 4)
        # load the ground truth for the bounding box
        truth = NP.zeros(number, dtype=_shortstruct)
        truth[0] = (78,  177, 177, 79)
        truth[1] = (334, 433, 177, 79)
        truth[2] = (78,  177, 433, 335)
        truth[3] = (334, 433, 433, 335)
        match = (truth==measures).all()
        assert_equal(match, True)
        # load the ground truth for the bounding box test image mean value
        voxel_truth = NP.zeros(number, dtype=NP.float64)
        voxel_truth = (80.0, 90.0, 100.0, 110.0)
        match = (voxel_truth==voxel_means).all()
        assert_equal(match, True)

        return

    def test_texture1(self):
        # [1] texture1 is delta functions and confirm the
        #     filter result is outer product of the L kernel
        M, Laws_LL = run_texture1()
        match = (Laws_LL==M).all()
        assert_equal(match, True)
        return

    def test_texture2(self):
        # [2] texture2 is 2 plane waves and assert the 20-element feature
        #     vector for each disc is correct
        tem = run_texture2()
        tem0 = tem[0]
        tem1 = tem[1]
        truth_tem0 = NP.array(
                        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                          0.        ,  0.13306101,  0.08511007,  0.05084148,  0.07550675,
                          0.4334695 ,  0.03715914,  0.00289055,  0.02755581,  0.48142046,
                          0.03137803,  0.00671277,  0.51568902,  0.01795249,  0.49102375,  1.
                        ], dtype=NP.float32)
        truth_tem1 = NP.array(
                        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                          0.        ,  0.02970393,  0.00164266,  0.00922416,  0.01221788,
                          0.51485199,  0.03298925,  0.02212243,  0.01912871,  0.48350537,
                          0.01125561,  0.00826189,  0.49437219,  0.00526817,  0.49736592,  1.
                        ], dtype=NP.float32)

        assert_array_almost_equal(tem0, truth_tem0, decimal=6)
        assert_array_almost_equal(tem1, truth_tem1, decimal=6)

        return

if __name__ == "__main__":
    inittest.main()

