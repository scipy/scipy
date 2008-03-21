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

def run_texture():
    filter      = seg.build_2d_kernel(hiFilterCutoff=60.0)
    img         = seg.build_test_discs()
    disc_mask   = seg.pre_filter(img, filter, low_threshold=50, high_threshold=255,
		                 conv_binary=1)
    label_disc_mask, disc_mask_groups = seg.get_blobs(disc_mask)
    disc_ROI    = seg.get_blob_regions(label_disc_mask, disc_mask_groups)
    laws_kernel = seg.build_laws_kernel() 
    impulse     = seg.build_test_impulses()
    calib       = seg.texture_filter(impulse, label_disc_mask, laws_kernel,
		                     ROI=disc_ROI, verbose=1)
    kernels = calib[0]
    x = laws_kernel['coefficients'][0]
    m = NP.outer(x, x)
    m = m * 2
    laws_LL = kernels[0, 50-4:50+3, 50-3:50+4]
    return m, laws_LL 

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

    def test_texture(self):
        # generate 4 discs; two tests (two test images) 
	# [1] image 1 is delta functions and confirm the
	#     filter result is outer product of the L kernel
	# [2] image 2 is 4 plane waves and assert the 20-element feature
	#     vector for each disc is correct
    	M, Laws_LL = run_texture()
    	match = (Laws_LL==M).all()
    	assert_equal(match, True)
    	return

if __name__ == "__main__":
    inittest.main()

