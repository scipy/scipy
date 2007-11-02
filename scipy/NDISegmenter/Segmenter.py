import numpy as N
import NDI_Segmenter as S
import volumeInput as V
import struct

def ShenCastan(image, IIRFilter=0.8, scLow=0.3, window=7, lowThreshold=220+2048, highThreshold=600+2048, dust=16):
	datatype = [('L', 'i', 1), ('R', 'i', 1), ('T', 'i', 1), ('B', 'i', 1),
            	    ('Label', 'i', 1), ('Area', 'i', 1), ('cX', 'f', 1), ('cY', 'f', 1),
            	    ('curveClose', 'i', 1), ('cXB', 'f', 1), ('cYB', 'f', 1), ('bLength', 'f', 1),
            	    ('minRadius', 'f', 1), ('maxRadius', 'f', 1), ('aveRadius', 'f', 1), ('ratio', 'f', 1),
            	    ('compactness', 'f', 1), ('voxelMean', 'f', 1), ('voxelVar', 'f', 1), ('TEM', 'f', 20)]
	labeledEdges, numberObjects = S.ShenCastanEdges(scLow, IIRFilter, window, lowThreshold, highThreshold, image)
	# allocated struct array for edge object measures. for now just the rect bounding box
	ROIList = N.zeros(numberObjects, dtype=datatype)
	# return the bounding box for each connected edge
	S.SetObjectStats(labeledEdges, ROIList)
	return labeledEdges, ROIList[ROIList['Area']>dust]

def Sobel(image, sLow=0.3, tMode=1, lowThreshold=220+2048, highThreshold=600+2048, BPHigh=10.0, apearture=21, dust=16):
	datatype = [('L', 'i', 1), ('R', 'i', 1), ('T', 'i', 1), ('B', 'i', 1),
            	    ('Label', 'i', 1), ('Area', 'i', 1), ('cX', 'f', 1), ('cY', 'f', 1),
            	    ('curveClose', 'i', 1), ('cXB', 'f', 1), ('cYB', 'f', 1), ('bLength', 'f', 1),
            	    ('minRadius', 'f', 1), ('maxRadius', 'f', 1), ('aveRadius', 'f', 1), ('ratio', 'f', 1),
            	    ('compactness', 'f', 1), ('voxelMean', 'f', 1), ('voxelVar', 'f', 1), ('TEM', 'f', 20)]
	# get sobel edge points. return edges that are labeled (1..numberObjects)
	labeledEdges, numberObjects = S.SobelEdges(sLow, tMode, lowThreshold, highThreshold, BPHigh, apearture, image)
	# allocated struct array for edge object measures. for now just the rect bounding box
	ROIList = N.zeros(numberObjects, dtype=datatype)
	# return the bounding box for each connected edge
	S.GetObjectStats(labeledEdges, ROIList)
	# thin (medial axis transform) of the sobel edges as the sobel produces a 'band edge'
	S.MorphoThinFilt(labeledEdges, ROIList)
	return labeledEdges, ROIList[ROIList['Area']>dust]

def Canny(image, cSigma=1.0, cLow=0.5, cHigh=0.8, tMode=1, lowThreshold=220+2048, highThreshold=600+2048,
          BPHigh=10.0, apearture=21, dust=16):
	datatype = [('L', 'i', 1), ('R', 'i', 1), ('T', 'i', 1), ('B', 'i', 1),
            	    ('Label', 'i', 1), ('Area', 'i', 1), ('cX', 'f', 1), ('cY', 'f', 1),
            	    ('curveClose', 'i', 1), ('cXB', 'f', 1), ('cYB', 'f', 1), ('bLength', 'f', 1),
            	    ('minRadius', 'f', 1), ('maxRadius', 'f', 1), ('aveRadius', 'f', 1), ('ratio', 'f', 1),
            	    ('compactness', 'f', 1), ('voxelMean', 'f', 1), ('voxelVar', 'f', 1), ('TEM', 'f', 20)]
	# get canny edge points. return edges that are labeled (1..numberObjects)
	labeledEdges, numberObjects = S.CannyEdges(cSigma, cLow, cHigh, tMode, lowThreshold, highThreshold, 
			                           BPHigh, apearture, image)
	# allocated struct array for edge object measures. for now just the rect bounding box
	ROIList = N.zeros(numberObjects, dtype=datatype)
	# return the bounding box for each connected edge
	S.GetObjectStats(labeledEdges, ROIList)
	return labeledEdges, ROIList[ROIList['Area']>dust]

def GetShapeMask(labeledEdges, ROIList):
	# pass in Sobel morph-thinned labeled edge image (LEI) and ROIList
	# GetShapeMask will augment the ROI list
	# labeledEdges is the original edge image and overwritten as mask image
	# maskImage is the mask that is used for blob texture / pixel features
	S.BuildBoundary(labeledEdges, ROIList)
	return 

def GetVoxelMeasures(rawImage, labeledEdges, ROIList):
	#
	# pass raw image, labeled mask and the partially filled ROIList
	# VoxelMeasures will fill the voxel features in the list
	#
	S.VoxelMeasures(rawImage, labeledEdges, ROIList)
	return 

def GetTextureMeasures(rawImage, labeledEdges, ROIList):
	#
	# pass raw image, labeled mask and the partially filled ROIList
	# VoxelMeasures will fill the texture (Law's, co-occurence, Gabor) features in the list
	#
	S.TextureMeasures(rawImage, labeledEdges, ROIList)
	return 

def SegmentRegions(volSlice=112):
	# get slice from the CT volume
    	image = GetSliceFromVolume(volSlice)
	# need a copy of original image as filtering will occur on the extracted slice
    	sourceImage = image.copy()
	# Sobel is the first level segmenter. Sobel magnitude and MAT (medial axis transform)
	# followed by connected component analysis. What is returned is labeled edges and the object list
    	labeledMask, ROIList = Sobel(image)
	# From the labeled edges and the object list get the labeled mask for each blob object
    	GetShapeMask(labeledMask, ROIList)
	# Use the labeled mask and source image (raw) to get voxel features 
    	GetVoxelMeasures(sourceImage, labeledMask, ROIList)
	# Use the labeled mask and source image (raw) to get texture features 
	GetTextureMeasures(sourceImage, labeledMask, ROIList)
	return sourceImage, labeledMask, ROIList

def GrowRegions(volSlice=112):
	# get slice from the CT volume
    	image = GetSliceFromVolume(volSlice)
	regionMask, numberRegions = RegionGrow(image)
	return regionMask, numberRegions 


def RegionGrow(image, lowThreshold=220+2048, highThreshold=600+2048, open=7, close=7):
	# morphology filters need to be clipped to 11 max and be odd
	regionMask, numberRegions = S.RegionGrow(lowThreshold, highThreshold, close, open, image)
	return regionMask, numberRegions
          

def GetSlice(imageName='junk.raw', bytes=2, rows=512, columns=512):
	# get a slice alrady extracted from the CT volume
	image = open(imageName, 'rb')
	slice = image.read(rows*columns*bytes)
	values = struct.unpack('h'*rows*columns, slice)
	ImageSlice = N.array(values, dtype=float).reshape(rows, columns)
	return (ImageSlice)

def GetSliceFromVolume(mySlice, rows=512, columns=512, bytes=2):
	# extract a slice the CT volume. Hardwirred 
	image = open('C:\PythonStuff\CardiacCT.vol', 'rb')
	image.seek(mySlice*rows*columns*bytes)
	slice = image.read(rows*columns*bytes)
	values = struct.unpack('h'*rows*columns, slice)
	ImageSlice = N.array(values, dtype=float).reshape(rows, columns)
	return (ImageSlice+2048)

def SaveSlice(mySlice, filename='junk.raw', rows=512, columns=512, bytes=2):
	# just save the slice to a fixed file
	slice = mySlice.astype(int)
	image = open(filename, 'wb')
	image.write(slice)
	image.close()
	return


