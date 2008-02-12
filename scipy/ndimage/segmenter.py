import math
import numpy as N
import scipy.ndimage._segment as S

# WARNING: _objstruct data structure mirrors a corresponding data structure
# in ndImage_Segmenter_structs.h that is built into the _segment.so library.
# These structs must match!  
_objstruct = N.dtype([('L', 'i'),
                    ('R', 'i'),
                    ('T', 'i'),
                    ('B', 'i'),
                    ('Label', 'i'),
                    ('Area', 'i'),
                    ('cX', 'f'),
                    ('cY', 'f'),
                    ('curveClose', 'i'),
                    ('cXB', 'f'),
                    ('cYB', 'f'),
                    ('bLength', 'f'),
                    ('minRadius', 'f'),
                    ('maxRadius', 'f'),
                    ('aveRadius', 'f'),
                    ('ratio', 'f'),
                    ('compactness', 'f'),
                    ('voxelMean', 'f'),
                    ('voxelVar', 'f'),
                    ('TEM', 'f', 20)]
                   )


def shen_castan(image, IIRFilter=0.8, scLow=0.3, window=7, lowThreshold=220+2048,
                highThreshold=600+2048, dust=16):
    """
        labeledEdges, ROIList = shen_castan(image, [default])

        implements Shen-Castan edge finding

        Inputs - image, IIR filter, shen_castan_low, window, low_threshold, high_threshold, dust 
        - image is the numarray 2D image
        - IIR filter is filter parameter for exponential filter
        - shen_castan_low is edge threshold is range (0.0, 1.0]
        - window is search window for edge detection
        - low_ and high_ threshold are density values 
        - dust is blob filter. blob area (length x width of bounding box) under this
          size threshold are filtered (referred to as dust and blown away)

        Outputs - labeledEdges, ROIList[>dust] 
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList[>dust] is a blob feature list. Only values
          with bounding box area greater than dust threshold are returned

    """
    labeledEdges, numberObjects = S.shen_castan_edges(scLow, IIRFilter, window, 
                                                      lowThreshold, highThreshold, image)
    # allocated struct array for edge object measures. for now just the rect bounding box
    ROIList = N.zeros(numberObjects, dtype=_objstruct)
    # return the bounding box for each connected edge
    S.get_object_stats(labeledEdges, ROIList)
    return labeledEdges, ROIList[ROIList['Area']>dust]

def sobel(image, sLow=0.3, tMode=1, lowThreshold=220+2048, highThreshold=600+2048, BPHigh=10.0, 
          apearture=21, dust=16):
    """
        labeledEdges, ROIList = sobel(image, [default])

        implements sobel magnitude edge finding

        Inputs - image, sobel_low, tMode, low_threshold, high_threshold, 
                        high_filter_cutoff, filter_aperature, dust
        - image is the numarray 2D image
        - sobel_low is edge threshold is range (0.0, 1.0]
        - tMode is threshold mode: 1 for ave, 2 for mode (histogram peak)
        - low_ and high_ threshold are density values 
        - high_filter_cutoff is digital high frequency cutoff in range (0.0, 180.0]
        - aperature is odd filter kernel length
        - dust is blob filter. blob area (length x width of bounding box) under this
          size threshold are filtered (referred to as dust and blown away)

        Outputs - labeledEdges, ROIList[>dust] 
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList[>dust] is a blob feature list. Only values
          with bounding box area greater than dust threshold are returned

    """
    # get sobel edge points. return edges that are labeled (1..numberObjects)
    labeledEdges, numberObjects = S.sobel_edges(sLow, tMode, lowThreshold, 
                                                highThreshold, BPHigh, apearture, image)
    # allocated struct array for edge object measures. for now just the rect bounding box
    ROIList = N.zeros(numberObjects, dtype=_objstruct)
    # return the bounding box for each connected edge
    S.get_object_stats(labeledEdges, ROIList)
    # thin (medial axis transform) of the sobel edges as the sobel produces a 'band edge'
    S.morpho_thin_filt(labeledEdges, ROIList)
    return labeledEdges, ROIList[ROIList['Area']>dust]

def canny(image, cSigma=1.0, cLow=0.5, cHigh=0.8, tMode=1, lowThreshold=220+2048, 
          highThreshold=600+2048, BPHigh=10.0, apearture=21, dust=16):
    """
        labeledEdges, ROIList = canny(image, [default])

        implements canny edge finding

        Inputs - image, DG_sigma, canny_low, canny_high, tMode, low_threshold,
                 high_threshold, high_filter_cutoff, filter_aperature, dust
        - image is the numarray 2D image
        - DG_sigma is Gaussain sigma for the derivative-of-gaussian filter
        - clow is low edge threshold is range (0.0, 1.0]
        - chigh is high edge threshold is range (0.0, 1.0]
        - tMode is threshold mode: 1 for ave, 2 for mode (histogram peak)
        - low_ and high_ threshold are density values 
        - high_filter_cutoff is digital high frequency cutoff in range (0.0, 180.0]
        - high_filter_cutoff is digital high frequency cutoff in range (0.0, 180.0]
        - aperature is odd filter kernel length
        - dust is blob filter. blob area (length x width of bounding box) under this
          size threshold are filtered (referred to as dust and blown away)

        Outputs - labeledEdges, ROIList[>dust] 
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList[>dust] is a blob feature list. Only values
          with bounding box area greater than dust threshold are returned

    """
    # get canny edge points. return edges that are labeled (1..numberObjects)
    labeledEdges, numberObjects = S.canny_edges(cSigma, cLow, cHigh, tMode, lowThreshold, highThreshold, 
                                               BPHigh, apearture, image)
    # allocated struct array for edge object measures. for now just the rect bounding box
    ROIList = N.zeros(numberObjects, dtype=_objstruct)
    # return the bounding box for each connected edge
    S.get_object_stats(labeledEdges, ROIList)
    return labeledEdges, ROIList[ROIList['Area']>dust]

def get_shape_mask(labeledEdges, ROIList):
    """
        get_shape_mask(labeledEdges, ROIList)

        takes labeled edge image plus ROIList (blob descriptors) and generates
        boundary shape features and builds labeled blob masks. 'labeledEdges' 
        is over-written by 'labeledMask'. Adds features to ROIList structure

        Inputs - labeledEdges, ROIList
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList is a blob feature list. 

        Output - no return. edge image input is over-written with mask image.
                            ROIList added to.

    """

    # pass in Sobel morph-thinned labeled edge image (LEI) and ROIList
    # GetShapeMask will augment the ROI list
    # labeledEdges is the original edge image and overwritten as mask image
    # maskImage is the mask that is used for blob texture / pixel features
    S.build_boundary(labeledEdges, ROIList)
    return 

def get_voxel_measures(rawImage, labeledEdges, ROIList):
    """
        get_voxel_measures(rawImage, labeledEdges, ROIList)

        takes raw 2D image, labeled blob mask and ROIList. computes voxel features
        (moments, histogram) for each blob. Adds features to ROIList structure.

        Inputs - rawImage, labeledEdges, ROIList
        - rawImage is the original source 2D image
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList is a blob feature list. 

        Output - no return. ROIList added to.

    """
    #
    # pass raw image, labeled mask and the partially filled ROIList
    # VoxelMeasures will fill the voxel features in the list
    #
    S.voxel_measures(rawImage, labeledEdges, ROIList)
    return 

def get_texture_measures(rawImage, labeledEdges, ROIList):
    """
        get_texture_measures(rawImage, labeledEdges, ROIList)

        takes raw 2D image, labeled blob mask and ROIList. computes 2D 
        texture features using 7x7 Law's texture filters applied 
        to segmented blobs. TEM (texture energy metric) is computed 
        for each Law's filter image and stored in TEM part of ROIList.

        Inputs - rawImage, labeledEdges, ROIList
        - rawImage is the original source 2D image
        - labeledEdges is boundary (edges) of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList is a blob feature list. 

        Output - no return. ROIList added to.
    """
    #
    # pass raw image, labeled mask and the partially filled ROIList
    # VoxelMeasures will fill the texture (Law's, sub-edges, co-occurence, Gabor) features in the list
    #
    S.texture_measures(rawImage, labeledEdges, ROIList)
    return 

def segment_regions(filename):
    """
        sourceImage, labeledMask, ROIList = segment_regions()

        Inputs - No Input

        Outputs - sourceImage, labeledMask, ROIList
        - sourceImage is raw 2D image (default cardiac CT slice for demo
        - labeledMask is mask of segmented 'blobs', 
          numerically labeled by blob number
        - ROIList is numerical Python structure of intensity, shape and 
          texture features for each blob

        High level script calls Python functions:
            get_slice()            - a cardiac CT slice demo file
            sobel()                - sobel magnitude edge finder,
                                     returns connected edges
            get_shape_mask()       - gets segmented blob boundary and mask 
                                     and shape features
            get_voxel_measures()   - uses masks get object voxel moment 
                                     and histogram features 
            get_texture_measures() - uses masks get object 2D texture features 
    """
    # get slice from the CT volume
    image = get_slice(filename)
    # need a copy of original image as filtering will occur on the extracted slice
    sourceImage = image.copy()
    # Sobel is the first level segmenter. Sobel magnitude and MAT (medial axis transform)
    # followed by connected component analysis. What is returned is labeled edges and the object list
    labeledMask, ROIList = sobel(image)
    # From the labeled edges and the object list get the labeled mask for each blob object
    get_shape_mask(labeledMask, ROIList)
    # Use the labeled mask and source image (raw) to get voxel features 
    get_voxel_measures(sourceImage, labeledMask, ROIList)
    # Use the labeled mask and source image (raw) to get texture features 
    get_texture_measures(sourceImage, labeledMask, ROIList)
    return sourceImage, labeledMask, ROIList

def grow_regions(filename):
    """
        regionMask, numberRegions = region_grow()
        Inputs - No Input
        Outputs - regionMask, numberRegions 
        - regionMask is the labeled segment masks from 2D image
        - numberRegions is the number of segmented blobs

        High level script calls Python functions:
            get_slice()      - a cardiac CT slice demo file
            region_grow()    - "grows" connected blobs. default threshold 
                                and morphological filter structuring element
    """
    # get slice from the CT volume
    image = get_slice(filename)
    regionMask, numberRegions = region_grow(image)
    return regionMask, numberRegions 


def region_grow(image, lowThreshold=220+2048, highThreshold=600+2048, open=7, close=7):
    """
        regionMask, numberRegions = region_grow(image, [defaults])

        Inputs - image, low_threshold, high_threshold, open, close
        - image is the numarray 2D image
        - low_ and high_ threshold are density values 
        - open is open morphology structuring element
          odd size. 0 to turn off. max is 11
        - close is close morphology structuring element
          odd size. 0 to turn off. max is 11

        Outputs - regionMask, numberRegions 
        - regionMask is the labeled segment masks from 2D image
        - numberRegions is the number of segmented blobs
    """
    # morphology filters need to be clipped to 11 max and be odd
    regionMask, numberRegions = S.region_grow(lowThreshold, highThreshold, close, open, image)
    return regionMask, numberRegions
      

def get_slice(imageName='slice112.raw', bytes=2, rows=512, columns=512):
    # get a slice alrady extracted from the CT volume
    #image = open(imageName, 'rb')
    #slice = image.read(rows*columns*bytes)
    #values = struct.unpack('h'*rows*columns, slice)
    #ImageSlice = N.array(values, dtype=float).reshape(rows, columns)

    ImageSlice = N.fromfile(imageName, dtype=N.uint16).reshape(rows, columns);

    # clip the ends for this test CT image file as the spine runs off the end of the image
    ImageSlice[505:512, :] = 0
    return (ImageSlice).astype(float)

def get_slice2(image_name='slice112.raw', bytes=2, shape=(512,512)):
    import mmap
    file = open(image_name, 'rb')
    mm = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
    slice = N.frombuffer(mm, dtype='u%d' % bytes).reshape(shape) 
    slice = slice.astype(float)
    # this is for the test CT as spine runs off back of image
    slice[505:512,:] = 0
    return slice

def save_slice(mySlice, filename='junk.raw', bytes=4):
    # just save the slice to a fixed file
    slice = mySlice.astype('u%d' % bytes)
    slice.tofile(filename)

def build_d_gauss_kernel(gWidth=21, sigma=1.0):

    """
    build the derivative of Gaussian kernel for Canny edge filter
    DGFilter = build_d_gauss_kernel(gWidth, sigma)
    Inputs:
        gWdith is width of derivative of Gaussian kernel
        sigma is sigma term of derivative of Gaussian kernel
    Output:
        DGFilter (a struct). Use in Canny filter call

    """
    kernel  = N.zeros((1+2*(gWidth-1)), dtype=float)
    indices = range(1, gWidth)  

    i = 0
    kernel[gWidth-1]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
    kernel[gWidth-1] *= -(i / (sigma * sigma))
    for i in indices:
        kernel[gWidth-1+i]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
        kernel[gWidth-1+i] *= -(i / (sigma * sigma))
        kernel[gWidth-1-i]  = -kernel[gWidth-1+i]

    DGFilter= {'kernelSize' : gWidth, 'coefficients': kernel} 

    return DGFilter

def build_2d_kernel(aperature=21, hiFilterCutoff=10.0):

    """
    build flat FIR filter with sinc kernel
    this is bandpass, but low cutoff is 0.0
    Use in Sobel and Canny filter edge find as image pre-process

    FIRFilter = build_2d_kernel(aperature, hiFilterCutoff)
    Inputs:
        aperature is number of FIR taps in sinc kernel 
        hiFilterCutoff is digital frequency cutoff in range (0.0, 180.0)
    Output:
        FIRFilter (a struct)

    """

    rad = math.pi / 180.0
    HalfFilterTaps = (aperature-1) / 2
    kernel = N.zeros((aperature), dtype=N.float32)
    LC = 0.0
    HC = hiFilterCutoff * rad 
    t2 = 2.0 * math.pi
    t1 = 2.0 * HalfFilterTaps + 1.0
    indices = range(-HalfFilterTaps, HalfFilterTaps+1, 1)  
    j = 0
    for i in indices:
        if i == 0:
            tLOW  = LC
            tHIGH = HC
        else:
            tLOW  = math.sin(i*LC)/i
            tHIGH = math.sin(i*HC)/i
        # Hamming window
        t3 = 0.54 + 0.46*(math.cos(i*t2/t1))
        t4 = t3*(tHIGH-tLOW)
        kernel[j] = t4
        j += 1

    # normalize the kernel
    sum = kernel.sum()
    kernel /= sum

    FIRFilter= {'kernelSize' : aperature, 'coefficients': kernel} 

    return FIRFilter


def build_laws_kernel():

    """
    build 6 length-7 Law's texture filter masks
    mask names are: 'L', 'S', 'E', 'W', 'R', 'O'

    LAWSFilter = build_laws_kernel()

    Inputs:
        None

    Output:
        LAWSFilter (a struct)

    """
    aperature = (6, 7)
    coefficients = N.zeros((aperature), dtype=N.float32)
    names = ('L', 'E', 'S', 'W', 'R', 'O' )

    coefficients[0, :] =  ( 1.0,  6.0,  15.0, 20.0,  15.0,  6.0,  1.0 )
    coefficients[1, :] =  (-1.0, -4.0,  -5.0,  0.0,   5.0,  4.0,  1.0 )
    coefficients[2, :] =  (-1.0, -2.0,   1.0,  4.0,   1.0, -2.0, -1.0 )
    coefficients[3, :] =  (-1.0,  0.0,   3.0,  0.0,  -3.0,  0.0,  1.0 )
    coefficients[4, :] =  ( 1.0, -2.0,  -1.0,  4.0,  -1.0, -2.0,  1.0 )
    coefficients[5, :] =  (-1.0,  6.0, -15.0, 20.0, -15.0,  6.0, -1.0 )

    LAWSFilter= {'numKernels' : 6, 'kernelSize' : 7, 'coefficients': coefficients, 'names': names} 

    return LAWSFilter

def build_morpho_thin_masks():

    """
    build 2 sets (J and K) of 8 3x3 morphology masks (structuring elements)
    to implement thinning (medial axis transformation - MAT)

    MATFilter = build_morpho_thin_masks()

    Inputs:
        None

    Output:
        MATFilter (a struct)

    """

    # (layers, rows, cols)
    shape  = (8, 3, 3)
    J_mask = N.zeros((shape), dtype=N.ushort)
    K_mask = N.zeros((shape), dtype=N.ushort)

    # load the 8 J masks for medial axis transformation
    J_mask[0][0][0] = 1;
    J_mask[0][0][1] = 1;
    J_mask[0][0][2] = 1;
    J_mask[0][1][1] = 1;

    J_mask[1][0][1] = 1;
    J_mask[1][1][1] = 1;
    J_mask[1][1][2] = 1;

    J_mask[2][0][0] = 1;
    J_mask[2][1][0] = 1;
    J_mask[2][2][0] = 1;
    J_mask[2][1][1] = 1;

    J_mask[3][0][1] = 1;
    J_mask[3][1][0] = 1;
    J_mask[3][1][1] = 1;

    J_mask[4][0][2] = 1;
    J_mask[4][1][1] = 1;
    J_mask[4][1][2] = 1;
    J_mask[4][2][2] = 1;

    J_mask[5][1][0] = 1;
    J_mask[5][1][1] = 1;
    J_mask[5][2][1] = 1;

    J_mask[6][1][1] = 1;
    J_mask[6][2][0] = 1;
    J_mask[6][2][1] = 1;
    J_mask[6][2][2] = 1;

    J_mask[7][1][1] = 1;
    J_mask[7][1][2] = 1;
    J_mask[7][2][1] = 1;


    # load the 8 K masks for medial axis transformation
    K_mask[0][2][0] = 1;
    K_mask[0][2][1] = 1;
    K_mask[0][2][2] = 1;

    K_mask[1][1][0] = 1;
    K_mask[1][2][0] = 1;
    K_mask[1][2][1] = 1;

    K_mask[2][0][2] = 1;
    K_mask[2][1][2] = 1;
    K_mask[2][2][2] = 1;

    K_mask[3][1][2] = 1;
    K_mask[3][2][1] = 1;
    K_mask[3][2][2] = 1;

    K_mask[4][0][0] = 1;
    K_mask[4][1][0] = 1;
    K_mask[4][2][0] = 1;

    K_mask[5][0][1] = 1;
    K_mask[5][0][2] = 1;
    K_mask[5][1][2] = 1;

    K_mask[6][0][0] = 1;
    K_mask[6][0][1] = 1;
    K_mask[6][0][2] = 1;

    K_mask[7][0][0] = 1;
    K_mask[7][0][1] = 1;
    K_mask[7][1][0] = 1;

    MATFilter = {'number3x3Masks' : 8, 'jmask' : J_mask, 'kmask' : K_mask} 

    return MATFilter 

