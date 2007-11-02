import volumeInput as V
import Segmenter as S
image = V.getFilteredSlice('slice112.raw')
sourceImage = image.copy()
edges, objects = S.Sobel(image)
S.GetShapeMask(edges, objects)
S.GetVoxelMeasures(sourceImage, edges, objects)
S.GetTextureMeasures(sourceImage, edges, objects)

import volumeInput as V
import Segmenter as S
sourceImage, labeledMask, ROIList = S.SegmentRegions()



image = S.GetSliceFromVolume()
edges, objects = S.Canny(image)

ShenCastanLow = 0.3
b             = 0.5
window        = 7
lowThreshold  = 220 + 2048
highThreshold = 600 + 2048
edges, groups = S.ShenCastan(image, ShenCastanLow, b, window, lowThreshold, highThreshold)


import Segmenter as S
image, mask, list = S.SegmentRegions()


import Segmenter as S
regionMask, numberRegions = S.GrowRegions()
regionMask.max()

S.SaveSlice(regionMask, 'regionMask.raw')

// for display of dumped .raw files using matplotlib


import volumeInput as V
rawslice = V.getFilteredSlice("source.raw", bytes=4)
edgeslice = V.getFilteredSlice("labeledMask.raw", bytes=4)

pylab.figure(1)
pylab.title('raw Image')
pylab.imshow(rawslice)
pylab.figure(3)
pylab.title('Edge Image')
pylab.imshow(edgeslice)


