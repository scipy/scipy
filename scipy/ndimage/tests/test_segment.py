
from scipy.testing import *
from scipy.ndimage._segmenter import *

inputname = 'slice112.raw'

from os.path import join, dirname
filename = join(dirname(__file__), inputname)

class TestSegment(TestCase):
    def test1(self):
        image = get_slice(filename)
        sourceImage = image.copy()
        edges, objects = sobel(image)
        get_shape_mask(edges, objects)
        get_voxel_measures(sourceImage, edges, objects)
        get_texture_measures(sourceImage, edges, objects)
	# measure the compactness and object boundry length
	# Ventricle measure
	assert_almost_equal(objects[7]['compactness'], 0.25657323, 4)
	assert_almost_equal(objects[7]['bLength'], 1215.70980000, 4)
	# Aorta measure
	assert_almost_equal(objects[13]['compactness'], 0.91137904, 4)
	assert_almost_equal(objects[13]['bLength'], 198.338090000, 4)

    def test2(self):
        sourceImage, labeledMask, ROIList = segment_regions(filename)
	# measure the compactness and object boundry length
	# Ventricle measure
	assert_almost_equal(ROIList[7]['compactness'], 0.25657323, 4)
	assert_almost_equal(ROIList[7]['bLength'], 1215.70980000, 4)
	# Aorta measure
	assert_almost_equal(ROIList[13]['compactness'], 0.91137904, 4)
	assert_almost_equal(ROIList[13]['bLength'], 198.338090000, 4)

    def test3(self):
        regionMask, numberRegions = grow_regions(filename)
	number_of_regions = regionMask.max()
	assert_equal(number_of_regions, 21)

    
if __name__ == "__main__":
    inittest.main()
