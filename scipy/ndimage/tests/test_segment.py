
from scipy.testing import *
from scipy.ndimage.segmenter import *

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

    def test2(self):
        sourceImage, labeledMask, ROIList = segment_regions(filename)

    def test3(self):
        regionMask, numberRegions = grow_regions(filename)
        regionMask.max()
        #save_slice(regionMask, 'regionMask.raw')

    
if __name__ == "__main__":
    inittest.main()
