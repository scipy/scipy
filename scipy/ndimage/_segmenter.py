import math
import numpy as N
import pylab as P
import scipy.ndimage._segment as S

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


# Issue warning regarding heavy development status of this module
import warnings
_msg = "The segmentation code is under heavy development and therefore the \
public API will change in the future.  The NIPY group is actively working on \
this code, and has every intention of generalizing this for the Scipy \
community.  Use this module minimally, if at all, until it this warning is \
removed."
warnings.warn(_msg, UserWarning)

def canny_hysteresis(magnitude, canny_stats):
    """
    edge_image = canny_hysteresis(magnitude, canny_stats)

    hystereis stage of Canny filter

    Parameters 
    ..........

    magnitude : {nd_array}
        the output from the canny_nonmax_supress() method

    canny_stats : {dictionary}
        contains the low and high thesholds determined from canny_nonmax_supress()

    Returns 
    ..........
    edge_image : {nd_array}
        the labeled edge image that can be displayed and used for later processing

    """
    [rows, cols] = magnitude.shape
    edge_image = N.zeros(rows*cols, dtype=N.int16).reshape(rows, cols)
    S.canny_hysteresis(magnitude, edge_image, canny_stats['low'], canny_stats['high'])

    return edge_image

def canny_nonmax_supress(horz_DGFilter, vert_DGFilter, img_means, thres=0.5, 
		         mode=1, canny_l=0.5, canny_h=0.8):
    """
    magnitude, canny_stats = canny_nonmax_supress(horz_DGFilter, vert_DGFilter, img_means,
		                          thres=0.5, mode=1, canny_l=0.5, canny_h=0.8)

    non-max supression stage of Canny filter

    Parameters 
    ..........

    horz_DGFilter : {nd_array}
        the horizonal filtered image using the derivative of Gaussian kernel filter.
	this is the output of the canny_filter method

    vert_DGFilter : {nd_array}
        the vertical filtered image using the derivative of Gaussian kernel filter.
	this is the output of the canny_filter method

    img_means : {dictionary}
        mean X and Y values of edge signals determined from canny_filter

    thres : {float}, optional
        edge threshold applied to x and y filtered means

    mode : {0, 1}, optional
        threshold based on histogram mean(0) or mode(1)

    canny_low : {float}, optional
        low threshold applied to magnitude filtered 
    
    canny_high : {float}, optional
        high threshold applied to magnitude filtered 

    Returns 
    ..........

    magnitude : {nd_array}
        magnitude of X and Y filtered for critical samples

    canny_stats : {dictionary}
        mean, low and high to be used for hysteresis

    """
    [rows, cols] = horz_DGFilter.shape
    magnitude = N.zeros(rows*cols, dtype=N.float64).reshape(rows, cols)
    aveMag, canny_low, canny_high = S.canny_nonmax_supress(horz_DGFilter, vert_DGFilter,
		                        magnitude, mode, img_means['x-dg']*thres,
				        img_means['y-dg']*thres, canny_l, canny_h)

    canny_stats = {'mean' : aveMag, 'low' : canny_low, 'high' : canny_high} 

    return magnitude, canny_stats

def canny_filter(slice, dg_kernel):
    """
    horz_DGFilter, vert_DGFilter, img_means = canny_filter(slice, dg_kernel)

    Canny derivative of Gaussian filter 
    returns the X and Y filterd image

    Parameters 
    ..........

    slice : {nd_array}
        2D image array

    dg_kernel : {dictionary}
        derivative of Gaussian kernel from build_d_gauss_kernel()

    Returns 
    ..........

    horz_DGFilter : {nd_array}
        X filtered image 

    vert_DGFilter : {nd_array}
        Y filtered image 

    img_means : {dictionary}


    """
    [rows, cols] = slice.shape
    horz_DGFilter = N.zeros(rows*cols, dtype=N.float64).reshape(rows, cols)
    vert_DGFilter = N.zeros(rows*cols, dtype=N.float64).reshape(rows, cols)
    aveX, aveY = S.canny_filter(slice, horz_DGFilter, vert_DGFilter,
		   dg_kernel['coefficients'], dg_kernel['kernelSize'])

    img_means = {'x-dg' : aveX, 'y-dg' : aveY} 

    return horz_DGFilter, vert_DGFilter, img_means


def mat_filter(label_image, thin_kernel, ROI=None):
    """
    mat_image = mat_filter(label_image, thin_kernel, ROI=None)

    takes the ROI dictionary with the blob bounding boxes and thins each blob
    giving the medial axis. if ROI is null will create a single ROI with the
    bounding box set equal to the full image

    Parameters 
    ..........

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    thin_kernel : {dictionary}
        set of 8 'J' and 'K' 3x3 masks from build_morpho_thin_masks() method

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes

    Returns 
    ..........

    mat_image : {nd_array}
        thinned edge image

    """
    if ROI==None:
        ROIList = N.zeros(1, dtype=_objstruct)
	[rows, cols] = label_image.shape
        ROIList['L'] = 2
        ROIList['R'] = cols-3
        ROIList['B'] = 2
        ROIList['T'] = rows-3

    [rows, cols] = label_image.shape
    # destination image
    thin_edge_image = N.zeros(rows*cols, dtype=N.uint16).reshape(rows, cols)
    mat_image       = N.zeros(rows*cols, dtype=N.uint16).reshape(rows, cols)
    # scratch memory for thin 
    input           = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)
    cinput          = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)
    erosion         = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)
    dialation       = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)
    hmt             = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)
    copy            = N.zeros(rows*cols, dtype=N.uint8).reshape(rows, cols)

    number_regions = ROI.size
    indices = range(0, number_regions)
    inflate = 1
    for i in indices:
	left     = ROI[i]['L']-1
	right    = ROI[i]['R']+1
	bottom   = ROI[i]['B']-1
	top      = ROI[i]['T']+1
	Label    = ROI[i]['Label']
	if left < 0: 
	    left = 0
	if bottom < 0: 
	    bottom = 0
	if right > cols-1: 
	    right = cols-1
	if top > rows-1: 
	    top = rows-1

	roi_rows = top-bottom+2*inflate
	roi_cols = right-left+2*inflate
	rgrows   = top-bottom
	rgcols   = right-left
	# clear the memory
	input[0:roi_rows, 0:roi_cols] = 0
	# load the labeled region 
	input[inflate:inflate+rgrows, inflate:inflate+rgcols] \
	     [label_image[bottom:top, left:right]==Label] = 1 
	# thin this region
        S.thin_filter(thin_kernel['jmask'], thin_kernel['kmask'], thin_kernel['number3x3Masks'],
		      roi_rows, roi_cols, cols, input, cinput, erosion, dialation, hmt, copy);

	# accumulate the images (do not over-write). for overlapping regions
	input[inflate:rgrows+inflate,inflate:rgcols+inflate] \
	     [input[inflate:rgrows+inflate,inflate:rgcols+inflate]==1] = Label 
	thin_edge_image[bottom:top,left:right] = thin_edge_image[bottom:top,left:right] + \
	                                         input[inflate:rgrows+inflate,inflate:rgcols+inflate] 
	    

    # accumulate overlaps set back to binary at later date
    mat_image[:, :] = thin_edge_image[:, :]

    return mat_image


def get_blob_regions(labeled_image, groups, dust=16):
    """
    ROIList = get_blob_regions(labeled_image, groups, dust=16)

    get the bounding box for each labelled blob in the image
    allocate the dictionary structure that gets filled in with later
    stage processing to add blob features.

    Parameters 
    ..........

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    groups : {int}
        number of blobs in image determined by get_blobs() method

    Returns 
    ..........

    ROIList : {dictionary}
        structure that has the bounding box and area of each blob



    """
    ROIList = N.zeros(groups, dtype=_objstruct)
    # return the bounding box for each connected edge
    S.get_blob_regions(labeled_image, ROIList)

    return ROIList[ROIList['Area']>dust]


def get_blobs(binary_edge_image):
    """

    labeled_edge_image, groups = get_blobs(binary_edge_image)

    get the total number of blobs in a 2D image and convert the binary
    image to labelled regions

    Parameters 
    ..........

    binary_edge_image : {nd_array}
        an binary image

    Returns 
    ..........

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    groups : {int}
        number of blobs in image determined by get_blobs() method

    """
    [rows, cols] = binary_edge_image.shape
    labeled_edge_image = N.zeros(rows*cols, dtype=N.uint16).reshape(rows, cols)
    groups = S.get_blobs(binary_edge_image, labeled_edge_image)

    return labeled_edge_image, groups


def sobel_edges(sobel_edge_image, sobel_stats, mode=1, sobel_threshold=0.3):
    """
    sobel_edge = sobel_edges(sobel_edge_image, sobel_stats, mode=1, sobel_threshold=0.3)
    take sobel-filtered image and return binary edges

    Parameters 
    ..........

    sobel_edge_image : {nd_array}
        edge-filtered image from sobel_image() method

    sobel_stats : {dictionary}
        mean and nonzero min, max of sobel filtering

    mode : {0, 1}, optional
        threshold based on histogram mean(0) or mode(1)

    sobel_threshold : {float}, optional
        low threshold applied to edge filtered image for edge generation

    Returns 
    ..........

    sobel_edge : {nd_array}
        binary edge-image

    """
    [rows, cols] = sobel_edge_image.shape
    sobel_edge = N.zeros(rows*cols, dtype=N.uint16).reshape(rows, cols)
    S.sobel_edges(sobel_edge_image, sobel_edge, sobel_stats['ave_gt0'], sobel_stats['min_gt0'],
                  sobel_stats['max_gt0'], mode, sobel_threshold)

    return sobel_edge


def sobel_image(filtered_slice):
    """
    sobel_edge_image, sobel_stats = sobel_image(filtered_slice)

    take 2D raw or filtered image slice and get sobel-filtered image 

    Parameters 
    ..........

    filtered_slice : {nd_array}
        raw or pre-processed (filtered and thresholded) 2D image 

    Returns 
    ..........

    sobel_edge_image : {nd_array}
        edge-filtered image from sobel_image() method

    sobel_stats : {dictionary}
        mean and nonzero min, max of sobel filtering

    """
    [rows, cols] = filtered_slice.shape
    sobel_edge_image = N.zeros(rows*cols, dtype=N.float64).reshape(rows, cols)
    pAve, min_value, max_value = S.sobel_image(filtered_slice, sobel_edge_image)

    # replace this with numpy calls for the image stats. but C ext is faster
    #   S.sobel_image(filtered_slice, sobel_edge_image)
    #   pAve = sobel_edge_image[sobel_edge_image>0].mean()
    #   min_value = sobel_edge_image[sobel_edge_image>0].min()
    #   max_value = sobel_edge_image[sobel_edge_image>0].max()

    sobel_stats= {'ave_gt0' : pAve, 'min_gt0': min_value, 'max_gt0': max_value} 

    return sobel_edge_image, sobel_stats

def pre_filter(slice, filter, low_threshold=2048+220, high_threshold=600+2048, conv_binary=0):
    """
    take 2D image slice and filter and pre-filter and threshold prior to segmentation

    """

    # make sure the input is 16 bits. this is input to edge machine
    # so can handle raw and 8 bit scaled inputs
    slice = slice.astype(N.int16)
    [rows, cols] = slice.shape
    edge_image = N.zeros(rows*cols, dtype=N.float64).reshape(rows, cols)
    S.edge_prefilter(low_threshold, high_threshold, filter['kernelSize'], filter['kernel'],
		     slice, edge_image)

    if conv_binary == 1:
        edge_image[edge_image>0] = 1
        edge_image = edge_image.astype(N.uint16)

    return edge_image

def get_max_bounding_box(ROI):
    max_area = ROI[:]['Area'].max()
    indices = range(0, ROI.size)
    for i in indices:
        if ROI[i]['Area'] == max_area:
	    left   = ROI[i]['L']
	    right  = ROI[i]['R']
	    top    = ROI[i]['T']
	    bottom = ROI[i]['B']

    bounding_box = {'left' : left, 'right' : right, 'top' : top, 'bottom' : bottom} 

    return bounding_box 

def set_draw_bounding_box(bounding_box):
    x = N.zeros(5, dtype=N.uint16)
    y = N.zeros(5, dtype=N.uint16)

    x[0] = bounding_box['left']
    x[1] = bounding_box['right']
    x[2] = bounding_box['right']
    x[3] = bounding_box['left']
    x[4] = bounding_box['left']

    y[0] = bounding_box['bottom']
    y[1] = bounding_box['bottom']
    y[2] = bounding_box['top']
    y[3] = bounding_box['top']
    y[4] = bounding_box['bottom']

    return x, y

def get_all_bounding_boxes(ROI):
    number = ROI.size
    indices = range(0, ROI.size)
    _shortstruct = N.dtype([('left', 'i'),
                            ('right', 'i'),
                            ('top', 'i'),
                            ('bottom', 'i')])
    measures = N.zeros(number, dtype=_shortstruct)
    for i in indices:
	measures[i]['left']   = ROI[i]['L']
	measures[i]['right']  = ROI[i]['R']
	measures[i]['top']    = ROI[i]['T']
	measures[i]['bottom'] = ROI[i]['B']

    return measures

def draw_all_bounding_boxes(measures):
    number = measures.size
    indices = range(0, measures.size)
    for i in indices:
        x, y = set_draw_bounding_box(measures[i])
	P.plot(x, y)

def build_test_discs():

    radius = 50
    rows   = 512
    cols   = 512
    test_image = N.zeros(rows*cols, dtype=N.int16).reshape(rows, cols)
    y_indices = N.array(range(-radius, radius+1))
    center_x = rows / 4
    center_y = cols / 4

    for i in y_indices:
	x = math.sqrt(float(radius)**2 - float(i)**2)
	test_image[1*center_y+i, 1*center_x-x:1*center_x+x] = 100
	test_image[1*center_y+i, 3*center_x-x:3*center_x+x] = 100
	test_image[3*center_y+i, 1*center_x-x:1*center_x+x] = 100
	test_image[3*center_y+i, 3*center_x-x:3*center_x+x] = 100

    return test_image

def get_slice(imageName='slice112.raw', bytes=2, rows=512, columns=512):
    # clip the ends for this test CT image file as the spine runs off the end of the image
    ImageSlice = N.fromfile(imageName, dtype=N.uint16).reshape(rows, columns)
    ImageSlice[505:512, :] = 0
    return ImageSlice

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
    kernel = N.zeros((aperature), dtype=N.float64)
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

    FIRFilter= {'kernelSize' : HalfFilterTaps, 'kernel': kernel} 

    return FIRFilter

def build_d_gauss_kernel(gWidth=20, sigma=1.0):

    """
    build the derivative of Gaussian kernel for Canny edge filter
    DGFilter = build_d_gauss_kernel(gWidth, sigma)
    Inputs:
        gWdith is width of derivative of Gaussian kernel
        sigma is sigma term of derivative of Gaussian kernel
    Output:
        DGFilter (a dictionary). Use in Canny filter call

    """

    kernel  = N.zeros((1+gWidth), dtype=N.float64)
    indices = range(0, gWidth)  

    for i in indices:
        kernel[i]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
        kernel[i] *= -(i / (sigma * sigma))

    DGFilter= {'kernelSize' : gWidth, 'coefficients': kernel} 

    return DGFilter

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
    size = (8*3*3)
    J_mask = N.zeros(size, dtype=N.int16)
    K_mask = N.zeros(size, dtype=N.int16)

    maskCols = 3
    # load the 8 J masks for medial axis transformation
    Column = 0
    J_mask[0+maskCols*(Column+0)] = 1
    J_mask[0+maskCols*(Column+1)] = 1
    J_mask[0+maskCols*(Column+2)] = 1
    J_mask[1+maskCols*(Column+1)] = 1

    Column += 3
    J_mask[0+maskCols*(Column+1)] = 1
    J_mask[1+maskCols*(Column+1)] = 1
    J_mask[1+maskCols*(Column+2)] = 1

    Column += 3
    J_mask[0+maskCols*(Column+0)] = 1
    J_mask[1+maskCols*(Column+0)] = 1
    J_mask[2+maskCols*(Column+0)] = 1
    J_mask[1+maskCols*(Column+1)] = 1

    Column += 3
    J_mask[0+maskCols*(Column+1)] = 1
    J_mask[1+maskCols*(Column+0)] = 1
    J_mask[1+maskCols*(Column+1)] = 1

    Column += 3
    J_mask[0+maskCols*(Column+2)] = 1
    J_mask[1+maskCols*(Column+1)] = 1
    J_mask[1+maskCols*(Column+2)] = 1
    J_mask[2+maskCols*(Column+2)] = 1

    Column += 3
    J_mask[1+maskCols*(Column+0)] = 1
    J_mask[1+maskCols*(Column+1)] = 1
    J_mask[2+maskCols*(Column+1)] = 1

    Column += 3
    J_mask[1+maskCols*(Column+1)] = 1
    J_mask[2+maskCols*(Column+0)] = 1
    J_mask[2+maskCols*(Column+1)] = 1
    J_mask[2+maskCols*(Column+2)] = 1

    Column += 3
    J_mask[1+maskCols*(Column+1)] = 1
    J_mask[1+maskCols*(Column+2)] = 1
    J_mask[2+maskCols*(Column+1)] = 1

    # load the 8 K masks for medial axis transformation
    Column = 0
    K_mask[2+maskCols*(Column+0)] = 1
    K_mask[2+maskCols*(Column+1)] = 1
    K_mask[2+maskCols*(Column+2)] = 1
    
    Column += 3
    K_mask[1+maskCols*(Column+0)] = 1
    K_mask[2+maskCols*(Column+0)] = 1
    K_mask[2+maskCols*(Column+1)] = 1
    
    Column += 3
    K_mask[0+maskCols*(Column+2)] = 1
    K_mask[1+maskCols*(Column+2)] = 1
    K_mask[2+maskCols*(Column+2)] = 1

    Column += 3
    K_mask[1+maskCols*(Column+2)] = 1
    K_mask[2+maskCols*(Column+1)] = 1
    K_mask[2+maskCols*(Column+2)] = 1

    Column += 3
    K_mask[0+maskCols*(Column+0)] = 1
    K_mask[1+maskCols*(Column+0)] = 1
    K_mask[2+maskCols*(Column+0)] = 1

    Column += 3
    K_mask[0+maskCols*(Column+1)] = 1
    K_mask[0+maskCols*(Column+2)] = 1
    K_mask[1+maskCols*(Column+2)] = 1

    Column += 3
    K_mask[0+maskCols*(Column+0)] = 1
    K_mask[0+maskCols*(Column+1)] = 1
    K_mask[0+maskCols*(Column+2)] = 1

    Column += 3
    K_mask[0+maskCols*(Column+0)] = 1
    K_mask[0+maskCols*(Column+1)] = 1
    K_mask[1+maskCols*(Column+0)] = 1

    MATFilter = {'number3x3Masks' : 8, 'jmask' : J_mask, 'kmask' : K_mask} 

    return MATFilter 

