import math
import numpy as NP
import scipy.ndimage._segment as S

_objstruct = NP.dtype([('Left', 'i'),
                       ('Right', 'i'),
                       ('Top', 'i'),
                       ('Bottom', 'i'),
                       ('Front', 'i'),
                       ('Back', 'i'),
                       ('Label', 'i'),
                       ('Mass', 'i'),
                       ('cX', 'f'),
                       ('cY', 'f'),
                       ('cZ', 'f'),
                       ('voxelMean', 'f'),
                       ('voxelVar', 'f'),
                       ('COM', 'f', 6),
                       ('TEM', 'f', 21)]
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
    ----------

    magnitude : {nd_array}
        the output from the canny_nonmax_supress() method

    canny_stats : {dictionary}
        contains the low and high thesholds determined from canny_nonmax_supress()

    Returns 
    ----------
    edge_image : {nd_array}
        the labeled edge image that can be displayed and used for later processing

    """
    [rows, cols] = magnitude.shape
    edge_image = NP.zeros(rows*cols, dtype=NP.int16).reshape(rows, cols)
    S.canny_hysteresis(magnitude, edge_image, canny_stats['low'], canny_stats['high'])

    return edge_image

def canny_nonmax_supress(horz_DGFilter, vert_DGFilter, img_means, thres=0.5, 
		         mode=1, canny_l=0.5, canny_h=0.8):
    """
    magnitude, canny_stats = canny_nonmax_supress(horz_DGFilter, vert_DGFilter, img_means,
		                          thres=0.5, mode=1, canny_l=0.5, canny_h=0.8)

    non-max supression stage of Canny filter

    Parameters 
    ----------

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
    ----------

    magnitude : {nd_array}
        magnitude of X and Y filtered for critical samples

    canny_stats : {dictionary}
        mean, low and high to be used for hysteresis

    """
    [rows, cols] = horz_DGFilter.shape
    magnitude = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
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
    ----------

    slice : {nd_array}
        2D image array

    dg_kernel : {dictionary}
        derivative of Gaussian kernel from build_d_gauss_kernel()

    Returns 
    ----------

    horz_DGFilter : {nd_array}
        X filtered image 

    vert_DGFilter : {nd_array}
        Y filtered image 

    img_means : {dictionary}


    """
    slice = slice.astype(NP.float64)
    [rows, cols] = slice.shape
    horz_DGFilter = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    vert_DGFilter = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    aveX, aveY = S.canny_filter(slice, horz_DGFilter, vert_DGFilter,
		   dg_kernel['coefficients'], dg_kernel['kernelSize'])

    img_means = {'x-dg' : aveX, 'y-dg' : aveY} 

    return horz_DGFilter, vert_DGFilter, img_means


def binary_edge(label_image, ROI):
    """
    binary_edge_image = binary_edge(label_image, ROI)

    takes the ROI dictionary with the blob bounding boxes and generates 
    the binary edge for each blob. The ROI mask is used to isolate
    the binary edges for later use (e.g. contour tracing).

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes

    Returns 
    ----------

    binary_edge_image : {nd_array}
        edge image for each ROI combined into a single image

    """

    [rows, cols]      = label_image.shape
    binary_edge_image = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    number_regions    = ROI.size
    indices           = range(0, number_regions)
    for i in indices:
	left   = ROI[i]['Left']-2
	right  = ROI[i]['Right']+2
	bottom = ROI[i]['Bottom']-2
	top    = ROI[i]['Top']+2
	Label  = ROI[i]['Label']
	if left < 0: 
	    left = 0
	if bottom < 0: 
	    bottom = 0
	if right > cols-1: 
	    right = cols-1
	if top > rows-1: 
	    top = rows-1

	roi_rows = top-bottom
	roi_cols = right-left
        label_region  = NP.zeros(roi_rows*roi_cols, dtype=NP.uint16).reshape(roi_rows, roi_cols)
        input = NP.zeros(roi_rows*roi_cols, dtype=NP.uint16).reshape(roi_rows, roi_cols)
	# load the labeled region 
	label_region[0:roi_rows, 0:roi_cols][label_image[bottom:top, left:right]==Label] = 1 
	S.binary_edge(label_region, input)
	input[0:roi_rows,0:roi_cols][input[0:roi_rows,0:roi_cols]==1] = Label
	binary_edge_image[bottom:top,left:right] = binary_edge_image[bottom:top,left:right] + \
	                                                       input[0:roi_rows,0:roi_cols] 

    return binary_edge_image


def roi_co_occurence(label_image, raw_image, ROI, distance=2, orientation=90, verbose=0):
    """
    roi_co_occurence(label_image, raw_image, ROI, distance=2, orientation=90, verbose=0)

    - OR -

    texture_arrays = roi_co_occurence(label_image, raw_image, ROI, distance=2, verbose=1)

    (N-S, E-W, NW-SE, NE-SW) computes one of 4 directional co-occurence matrices and features.
    In debug=1 will return the 4 joint histograms for each ROI. This is used to compute texture
    features for a pre-segmented region. The seg_co_occurence() method will be used for texture-
    based segmentation.

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    raw_image : {nd_array}
        raw image from which texture features get extracted 

    distance : {int}
        integer value of pixel offset in forming joint histogram. default value 2

    orientation : {45, 90, 135, 180}
        direction for pixel offet.

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes. The largest
	2D target bounding box is extracted.


    Returns 
    ----------

    co_occurence_images : {dictionary}
        contains 4 joint histogram images for each ROI 
	returned if verbose=1

    """

    if orientation != 45 and orientation != 90 and orientation != 135 and orientation != 180: 
        orientation = 90

    epsilon  = 2.2e-16 
    num_bits = 256
    copy_image = raw_image.copy()
    number_regions = ROI.size
    indices = range(0, number_regions)
    co_occurence_image_list = {}
    for i in indices:
        left   = ROI[i]['Left']
        right  = ROI[i]['Right']
        bottom = ROI[i]['Bottom']
        top    = ROI[i]['Top']
        Label  = ROI[i]['Label']
        rows   = top-bottom
        cols   = right-left
	# copy the mask to section image
        section = NP.zeros(rows*cols, dtype=label_image.dtype).reshape(rows, cols)
        section[0:rows, 0:cols][label_image[bottom:top, left:right]==Label] = 1
        source_region = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
        cocm_block = NP.zeros(num_bits*num_bits, dtype=NP.int32).reshape(num_bits, num_bits)
	source_region[0:rows, 0:cols] = copy_image[bottom:top, left:right] 
        # scale segment to 8 bits. this needs to be smarter (e.g. use integrated histogram method)
        max_value = source_region.max()
        min_value = source_region.min()
        scale = 255.0 / (max_value-min_value)
        image_roi = (scale*(source_region-min_value)).astype(NP.int16)
	# image_roi is short type
	S.roi_co_occurence(section, image_roi, cocm_block, distance, orientation)
        co_occurence_image_list[i] = cocm_block
	# normalize the joint histogram prior to feature extraction
	joint_histogram  = cocm_block.astype(NP.float64) 
	joint_histogram  = joint_histogram / joint_histogram.sum()
	# to prevent log(0)
	joint_histogram += epsilon
	# compute the com features
	energy = joint_histogram.std()
	H = joint_histogram * NP.log(joint_histogram)
	entropy = H.sum()
	r, c = joint_histogram.shape
	[a, b] = NP.mgrid[1:c+1, 1:r+1]
	contrast = ((NP.square(a-b))*joint_histogram).sum()
	d = 1.0 + NP.abs(a-b)
	homogeneity = (joint_histogram / d).sum()
	ROI[i]['COM'][0] = distance
	ROI[i]['COM'][1] = orientation 
	ROI[i]['COM'][2] = energy
	ROI[i]['COM'][3] = entropy
	ROI[i]['COM'][4] = contrast
	ROI[i]['COM'][5] = homogeneity

    if verbose == 1:
        return co_occurence_image_list
    else:
	return



def region_grow(label_image, raw_image, ROI, roi_index, roi_inflate,
		low_thresh=0.5, high_thresh=1.5, N_connectivity=3, debug=0):
    """
    region_grow(label_image, raw_image, ROI, roi_index, roi_inflate, stop_thresh)

    starting from a seed (masks or regions in the label_image) grow the regions based
    on (1) connectivity (2D or 3D) and (2) raw voxel values > stop threshold criterion.

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    raw_image : {nd_array}
        raw image from which texture features get extracted 

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes. The largest
	2D target bounding box is extracted.

    roi_index : {int}
        the single ROI element to apply region growing to.

    roi_inflate : {list}
        the maximum increase in the ROI bounding box. For 3D the tuple is [layers, rows, cols]
	and for 2D it is [rows, cols].

    low_thresh : {float}
        this is the percent of the voxel mean that the growing region must be GREATER than.
	region growing terminates when the raw_image is BELOW this value.

    high_thresh : {float}
        this is the percent of the voxel mean that the growing region must be LESS than.
	region growing terminates when the raw_image is ABOVE this value.

    N_connectivity : {int}
        for growing this indicates how connected in a 3x3 or 3x3x3 window the un-labeled
	sample is. Make less than full connected for growing boundaries

    Returns 
    ----------

    label : {nd_array}
        the label image with the selected ROI after region growing. only returned
	in debug mode.

    """

    _c_ext_struct = NP.dtype([('Left', 'i'),
                              ('Right', 'i'),
                              ('Top', 'i'),
                              ('Bottom', 'i'),
                              ('Front', 'i'),
                              ('Back', 'i'),
                              ('Label', 'i'),
                              ('Mass', 'i'),
                              ('cX', 'f'),
                              ('cY', 'f'),
                              ('cZ', 'f')]
                             )

    expanded_ROI = NP.zeros(1, dtype=_c_ext_struct)

    dimensions  = label_image.ndim

    if dimensions == 3:
        z_ext = roi_inflate[0]
        y_ext = roi_inflate[1]
        x_ext = roi_inflate[2]
        [layers, rows, cols]  = label_image.shape
    else:
        y_ext = roi_inflate[0]
        x_ext = roi_inflate[1]
        [rows, cols]  = label_image.shape

    if dimensions == 2:  
        left    = ROI[roi_index]['Left']-x_ext
        right   = ROI[roi_index]['Right']+x_ext
        bottom  = ROI[roi_index]['Bottom']-y_ext
        top     = ROI[roi_index]['Top']+y_ext
        Label   = ROI[roi_index]['Label']
        lcutoff = low_thresh  * ROI[roi_index]['voxelMean']
        hcutoff = high_thresh * ROI[roi_index]['voxelMean']
   	if left < 0: 
           left = 0
    	if bottom < 0: 
            bottom = 0
    	if right > cols-1: 
            right = cols-1
    	if top > rows-1: 
            top = rows-1
        expanded_ROI['Left']   = left 
        expanded_ROI['Right']  = right 
        expanded_ROI['Top']    = top 
        expanded_ROI['Bottom'] = bottom 
        expanded_ROI['Label']  = Label 
	rows    = top-bottom
	cols    = right-left
        label   = NP.zeros(rows*cols, dtype=NP.int16).reshape(rows, cols)
        section = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
	label   = label_image[bottom:top, left:right].copy()
	section = (raw_image[bottom:top, left:right].astype(NP.float64)).copy()
    elif dimensions == 3:  
        left    = ROI[roi_index]['Left']-x_ext
        right   = ROI[roi_index]['Right']+x_ext
        bottom  = ROI[roi_index]['Bottom']-y_ext
        top     = ROI[roi_index]['Top']+y_ext
        front   = ROI[roi_index]['Front']-z_ext
        back    = ROI[roi_index]['Back']+z_ext
        Label   = ROI[roi_index]['Label']
        lcutoff = low_thresh  * ROI[roi_index]['voxelMean']
        hcutoff = high_thresh * ROI[roi_index]['voxelMean']
    	if left < 0: 
            left = 0
    	if bottom < 0: 
            bottom = 0
    	if right > cols-1: 
            right = cols-1
    	if top > rows-1: 
            top = rows-1
    	if front < 0: 
            front = 0
    	if back > layers-1: 
            back = layers-1
        expanded_ROI['Left']   = left 
        expanded_ROI['Right']  = right 
        expanded_ROI['Top']    = top 
        expanded_ROI['Bottom'] = bottom 
        expanded_ROI['Back']   = back 
        expanded_ROI['Front']  = front 
        expanded_ROI['Label']  = Label 
	rows    = top-bottom
	cols    = right-left
	layers  = back-front
        label   = NP.zeros(layers*rows*cols, dtype=NP.int16).reshape(layers, rows, cols)
	label   = label_image[front:back, bottom:top, left:right].copy()
        section = NP.zeros(layers*rows*cols, dtype=NP.float64).reshape(layers, rows, cols)
	section = (raw_image[front:back, bottom:top, left:right].astype(NP.float64)).copy()

    #
    # this newgrow_ROI gets filled in and the label image is grown
    #

    newgrow_ROI = NP.zeros(1, dtype=_c_ext_struct)
    S.region_grow(section, label, expanded_ROI, newgrow_ROI, lcutoff, hcutoff, Label, N_connectivity)

    if debug==1:  
	#
	# do not update ROI for index and the label_image 
	#
        return label

    else:
	#
	# update (overwrite) ROI for index and the label_image 
	#
        if dimensions == 2:  
	    ROI[roi_index]['Left']   = newgrow_ROI['Left']
	    ROI[roi_index]['Right']  = newgrow_ROI['Right']
	    ROI[roi_index]['Top']    = newgrow_ROI['Top']
	    ROI[roi_index]['Bottom'] = newgrow_ROI['Bottom']
	    left   = ROI[roi_index]['Left']
	    right  = ROI[roi_index]['Right']
	    top    = ROI[roi_index]['Top']
	    bottom = ROI[roi_index]['Bottom']
	    rows   = top-bottom
	    cols   = right-left
	    label_image[bottom:top,left:right] = label[0:rows,0:cols]
        elif dimensions == 3:  
            ROI[roi_index]['Left']   = newgrow_ROI['Left']
            ROI[roi_index]['Right']  = newgrow_ROI['Right']
            ROI[roi_index]['Top']    = newgrow_ROI['Top']
            ROI[roi_index]['Bottom'] = newgrow_ROI['Bottom']
            ROI[roi_index]['Front']  = newgrow_ROI['Front']
            ROI[roi_index]['Back']   = newgrow_ROI['Back']
	    left   = expanded_ROI['Left']
	    right  = expanded_ROI['Right']
	    top    = expanded_ROI['Top']
	    bottom = expanded_ROI['Bottom']
	    front  = expanded_ROI['Front']
	    back   = expanded_ROI['Back']
	    rows   = top-bottom
	    cols   = right-left
	    layers = back-front
	    label_image[front:back,bottom:top,left:right] = label[0:layers,0:rows,0:cols]
			     
        return 



def seg_co_occurence(raw_image, window=16, distance=2, orientation=90):
    """
    cocm_images = seg_co_occurence(raw_image, window=16, distance=2, orientation=90)

    (N-S, E-W, NW-SE, NE-SW) computes one of 4 directional co-occurence matrices and features.
    In debug=1 will return the 4 joint histograms for each ROI.

    The seg_co_occurence() method is used for texture-based segmentation. Feature images are
    returned from which segmentation can be later performed.

    ****
    NOTE: This is very slow and a fast method using Unsers histogram approximation will be 
    added in the future.
    ****


    Parameters 
    ----------

    raw_image : {nd_array}
        raw image from which texture features get extracted 

    window : {int}
        integer value of moving 2D window. Window slides in 2D over image and is the
	region-of-interest from which co-occurence texture features are extracted. The
	window is 2D square so only a single value is entered. Default window is 32x32. 

    distance : {int}
        integer value of pixel offset in forming joint histogram. default value 2

    orientation : {45, 90, 135, 180}
        direction for pixel offet.

    Returns 
    ----------

    cocm_images : {dictionary}
        
	co_occurence_feature_images. contains 4 normalized feature
	windows with keys: energy, entropy, contrast and homogeneity.

    """

    if orientation != 45 and orientation != 90 and orientation != 135 and orientation != 180: 
        orientation = 90

    epsilon      = 2.2e-16 
    num_bits     = 256
    copy_image   = raw_image.copy()
    [rows, cols] = copy_image.shape
    row_indices  = range(window, rows-window)
    col_indices  = range(window, cols-window)

    # create a fixed mask and scratch window for raw source
    section       = NP.ones(2*window*2*window, dtype=NP.int16).reshape(2*window, 2*window)
    source_region = NP.zeros(2*window*2*window, dtype=NP.float64).reshape(2*window, 2*window)

    # output images
    energy_image      = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    entropy_image     = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    homogeneity_image = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    contrast_image    = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    cocm_block        = NP.zeros(num_bits*num_bits, dtype=NP.int32).reshape(num_bits, num_bits)
    
    for i in row_indices:
        bottom = i - window
        top    = i + window
        for j in col_indices:
            left  = j - window
            right = j + window 
	    source_region[0:2*window, 0:2*window] = copy_image[bottom:top, left:right] 
            # scale segment to 8 bits. this needs to be smarter (e.g. use integrated histogram method)
            max_value = source_region.max()
            min_value = source_region.min()
            scale     = 255.0 / (max_value-min_value)
            image_roi = (scale*(source_region-min_value)).astype(NP.int16)
	    # image_roi is short type
	    cocm_block[:] = 0.0
	    S.roi_co_occurence(section, image_roi, cocm_block, distance, orientation)
	    # normalize the joint histogram prior to feature extraction
	    joint_histogram = cocm_block.astype(NP.float64) 
	    joint_histogram = joint_histogram / joint_histogram.sum()
	    # to prevent log(0)
	    joint_histogram += epsilon
	    # compute the com features
	    energy      = joint_histogram.std()
	    H           = joint_histogram * NP.log(joint_histogram)
	    entropy     = H.sum()
	    r, c        = joint_histogram.shape
	    [a, b]      = NP.mgrid[1:c+1, 1:r+1]
	    contrast    = ((NP.square(a-b))*joint_histogram).sum()
	    d           = 1.0 + NP.abs(a-b)
	    homogeneity = (joint_histogram / d).sum()
	    # store the feature pixel for the 4 images
	    energy_image[i, j]      = energy
	    entropy_image[i, j]     = entropy
	    contrast_image[i, j]    = contrast
	    homogeneity_image[i, j] = homogeneity

    scale_energy      = 1.0 / max(energy_image.max(), abs(energy_image.min()))
    scale_entropy     = 1.0 / max(entropy_image.max(), abs(entropy_image.min()))
    scale_contrast    = 1.0 / max(contrast_image.max(), abs(contrast_image.min()))
    scale_homogeneity = 1.0 / max(homogeneity_image.max(), abs(homogeneity_image.min()))

    energy_image      = scale_energy      * energy_image
    entropy_image     = scale_entropy     * entropy_image
    homogeneity_image = scale_homogeneity * homogeneity_image
    contrast_image    = scale_contrast    * contrast_image

    cocm_images = {'energy_image' : energy_image,  'entropy_image' : entropy_image, 
                   'homogeneity_image' : homogeneity_image,  'contrast_image' : contrast_image} 

    return cocm_images 


def roi_mat_filter(label_image, thin_kernel, ROI):
    """
    thin_edge_image = roi_mat_filter(label_image, thin_kernel, ROI)

    gets the largest object in the ROI list and returns the medial axis for
    that object. Idea is that the largest object is a reference, e.g. the skull
    for anatomical MRI.

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    thin_kernel : {dictionary}
        set of 8 'J' and 'K' 3x3 masks from build_morpho_thin_masks() method

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes. The largest
	2D target bounding box is extracted.

    Returns 
    ----------

    thin_edge_image : {nd_array}
        thinned edge image for the largest object.

    """


    [rows, cols] = label_image.shape
    # destination image
    thin_edge_image = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    # scratch memory for thin 
    input           = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    cinput          = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    erosion         = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    dialation       = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    hmt             = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    copy            = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)

    bbox = get_max_bounding_box(ROI)

    left     = bbox['Left']-1
    right    = bbox['Right']+1
    bottom   = bbox['Bottom']-1
    top      = bbox['Top']+1
    Label    = bbox['Label']

    if left < 0: 
        left = 0
    if bottom < 0: 
        bottom = 0
    if right > cols-1: 
        right = cols-1
    if top > rows-1: 
        top = rows-1

    inflate  = 1
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
	          roi_rows, roi_cols, cols, input, cinput, erosion, dialation, hmt, copy)

    # accumulate the images (do not over-write). for overlapping regions
    input[inflate:rgrows+inflate,inflate:rgcols+inflate] \
         [input[inflate:rgrows+inflate,inflate:rgcols+inflate]==1] = Label 
    thin_edge_image[bottom:top,left:right] = input[inflate:rgrows+inflate,inflate:rgcols+inflate] 

    return thin_edge_image

def mat_filter(label_image, thin_kernel, ROI=None):
    """
    mat_image = mat_filter(label_image, thin_kernel, ROI=None)

    takes the ROI dictionary with the blob bounding boxes and thins each blob
    giving the medial axis. if ROI is null will create a single ROI with the
    bounding box set equal to the full image

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    thin_kernel : {dictionary}
        set of 8 'J' and 'K' 3x3 masks from build_morpho_thin_masks() method

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes

    Returns 
    ----------

    mat_image : {nd_array}
        thinned edge image

    """
    if ROI==None:
        ROIList = NP.zeros(1, dtype=_objstruct)
	[rows, cols] = label_image.shape
        ROIList['Left']   = 2
        ROIList['Right']  = cols-3
        ROIList['Bottom'] = 2
        ROIList['Top']    = rows-3

    [rows, cols] = label_image.shape
    # destination image
    thin_edge_image = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    mat_image       = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    # scratch memory for thin 
    input           = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    cinput          = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    erosion         = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    dialation       = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    hmt             = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)
    copy            = NP.zeros(rows*cols, dtype=NP.uint8).reshape(rows, cols)

    number_regions = ROI.size
    indices = range(0, number_regions)
    inflate = 1
    for i in indices:
	left     = ROI[i]['Left']-1
	right    = ROI[i]['Right']+1
	bottom   = ROI[i]['Bottom']-1
	top      = ROI[i]['Top']+1
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
		      roi_rows, roi_cols, cols, input, cinput, erosion, dialation, hmt, copy)

	# accumulate the images (do not over-write). for overlapping regions
	input[inflate:rgrows+inflate,inflate:rgcols+inflate] \
	     [input[inflate:rgrows+inflate,inflate:rgcols+inflate]==1] = Label 
	thin_edge_image[bottom:top,left:right] = thin_edge_image[bottom:top,left:right] + \
	                                         input[inflate:rgrows+inflate,inflate:rgcols+inflate] 
	    

    # accumulate overlaps set back to binary at later date
    mat_image[:, :] = thin_edge_image[:, :]

    return mat_image


def laws_texture_filter(raw_image, label_image, laws_kernel, ROI=None, dc_thres=1.0,
		        mean_feature=1, verbose=0):
    """
    texture_images = laws_texture_filter(raw_image, label_image, laws_kernel, ROI=None, verbose=1)
    .
        OR
    .
    laws_texture_filter(raw_image, label_image, laws_kernel, ROI=None, verbose=0)

    Parameters 
    ----------

    raw_image : {nd_array}
        raw double image 

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    laws_kernel : {dictionary}
        set of 6 length-7 Law's texture feature kernels 

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes

    dc_thres : {float}
        used as a filter. Sets texture feature to 0.0 when the 
	mean level is above this. Removes the low frequency, high amplitude
	image regions from the feature list

    mean_feature : {0, 1}, optional
        when set to 1, the feature is the mean value of the
	selected Law's texture filter. When 0 the feature is
	the standard deviation.

    verbose : {0, 1}, optional
        determines if return is to include Law's filter images

    Returns 
    ----------

    laws_image : {dictionary}
        contains 21 Laws filtered  regions for each ROI 
	returned if verbose=1
        

    """
    if ROI==None:
        ROI= NP.zeros(1, dtype=_objstruct)
	[rows, cols] = label_image.shape
        ROI['Left']   = 2
        ROI['Right']  = cols-3
        ROI['Bottom'] = 2
        ROI['Top']    = rows-3

    laws_image_list = {}
    number_regions  = ROI.size
    layers          = laws_kernel['filters']
    indices         = range(0, number_regions)
    filters         = range(0, layers)
    for i in indices:
	left   = ROI[i]['Left']
	right  = ROI[i]['Right']
	bottom = ROI[i]['Bottom']
	top    = ROI[i]['Top']
	Label  = ROI[i]['Label']
	rows   = top-bottom
	cols   = right-left
        label_region  = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
        source_region = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
        laws_block    = NP.zeros(layers*rows*cols, dtype=NP.float32).reshape(layers, rows, cols)
	# load the labeled region 
	label_region[0:rows,  0:cols][label_image[bottom:top, left:right]==Label] = 1 
	source_region[0:rows, 0:cols] = raw_image[bottom:top, left:right] 

	S.laws_texture_metric(label_region, source_region, laws_block, laws_kernel['numKernels'],
		              laws_kernel['kernelSize'], laws_kernel['filters'],
		              laws_kernel['coefficients'][0], laws_kernel['coefficients'][1],
		              laws_kernel['coefficients'][2], laws_kernel['coefficients'][3],
		              laws_kernel['coefficients'][4], laws_kernel['coefficients'][5])

        for j in filters:
	    # compute the energy measure for each filter in the ROI
	    mask_image = laws_block[j, :, :][label_region[:, :]>0]
	    mean = abs(mask_image.mean())
	    std  = mask_image.std()
	    if mean > dc_thres:
	        mean = 0.0
	        std = 0.0
	    if mean_feature == 1:
	        ROI[i]['TEM'][j] = mean 
            else:
	        ROI[i]['TEM'][j] = std 

	ROI[i]['TEM'][:] = ROI[i]['TEM'][:] / ROI[i]['TEM'][:].max() 
        # accumulate the 21 Law's filtered ROI's and optional
	# return as image (3D)
        laws_image_list[i] = laws_block

    if verbose == 1:
	return laws_image_list
    else:
        return 



def get_voxel_measures(label_image, raw_image, ROI=None):
    """
    mat_image = mat_filter(label_image, raw_image, ROI=None)

    takes the ROI dictionary with the blob bounding boxes and gets the voxel measures
    from each ROI in the raw data.

    Parameters 
    ----------

    label_image : {nd_array}
        an image with labeled regions from get_blobs() method

    raw_image : {nd_array}
        the original double image (raw voxels) from which voxel measures are made

    ROI : {dictionary}
        Region of Interest structure that has blob bounding boxes


    Returns 
    ----------

    none

    """

    dimensions = label_image.ndim

    if ROI==None:
        ROIList = NP.zeros(1, dtype=_objstruct)
        if dimensions == 2:  
	    [rows, cols] = label_image.shape
            ROIList['Left']   = 1
            ROIList['Right']  = cols-1
            ROIList['Bottom'] = 1
            ROIList['Top']    = rows-1
        elif dimensions == 3:  
	    [layers, rows, cols] = label_image.shape
            ROIList['Left']   = 1
            ROIList['Right']  = cols-1
            ROIList['Bottom'] = 1
            ROIList['Top']    = rows-1
            ROIList['Front']  = 1
            ROIList['Back']   = layers-1

    number_regions = ROI.size
    indices = range(0, number_regions)
    inflate = 1
    for i in indices:
        if dimensions == 2:  
	    left   = ROI[i]['Left']
	    right  = ROI[i]['Right']
	    bottom = ROI[i]['Bottom']
	    top    = ROI[i]['Top']
	    Label  = ROI[i]['Label']
	    rows   = top-bottom-1
	    cols   = right-left-1
            section= NP.zeros(rows*cols, dtype=raw_image.dtype).reshape(rows, cols)
	    section = raw_image[bottom:top, left:right] \
	                       [label_image[bottom:top, left:right]==Label]
        elif dimensions == 3:  
	    left   = ROI[i]['Left']
	    right  = ROI[i]['Right']
	    bottom = ROI[i]['Bottom']
	    top    = ROI[i]['Top']
	    front  = ROI[i]['Front']
	    back   = ROI[i]['Back']
	    Label  = ROI[i]['Label']
	    rows   = top-bottom-1
	    cols   = right-left-1
	    layers = back-front-1
            section= NP.zeros(layers*rows*cols, dtype=raw_image.dtype).reshape(layers, rows, cols)
	    section = raw_image[front:back, bottom:top, left:right] \
			       [label_image[front:back, bottom:top, left:right]==Label]

	mask = section[section>0]
	ROI[i]['voxelMean'] = mask.mean()
	ROI[i]['voxelVar']  = mask.std()

    return 


def get_blob_regions(labeled_image, groups, dust=16):
    """
    ROIList = get_blob_regions(labeled_image, groups, dust=16)

    get the bounding box for each labelled blob in the image
    allocate the dictionary structure that gets filled in with later
    stage processing to add blob features.

    Parameters 
    ----------

    label_image : {nd_array}
        a 2D or 3D image with labeled regions from get_blobs() method

    groups : {int}
        number of blobs in image determined by get_blobs() method

    Returns 
    ----------

    ROIList : {dictionary}
        structure that has the bounding box and area of each blob



    """

    _c_ext_struct = NP.dtype([('Left', 'i'),
                              ('Right', 'i'),
                              ('Top', 'i'),
                              ('Bottom', 'i'),
                              ('Front', 'i'),
                              ('Back', 'i'),
                              ('Label', 'i'),
                              ('Mass', 'i'),
                              ('cX', 'f'),
                              ('cY', 'f'),
                              ('cZ', 'f')]
                             )

    c_ext_ROI = NP.zeros(groups, dtype=_c_ext_struct)
    ROIList = NP.zeros(groups, dtype=_objstruct)
    # return the bounding box for each connected edge
    S.get_blob_regions(labeled_image, c_ext_ROI)

    indices = range(0, groups)
    for i in indices:
	ROIList[i]['Left']   = c_ext_ROI[i]['Left']
	ROIList[i]['Right']  = c_ext_ROI[i]['Right']
	ROIList[i]['Bottom'] = c_ext_ROI[i]['Bottom']
	ROIList[i]['Top']    = c_ext_ROI[i]['Top']
	ROIList[i]['Front']  = c_ext_ROI[i]['Front']
	ROIList[i]['Back']   = c_ext_ROI[i]['Back']
	ROIList[i]['Label']  = c_ext_ROI[i]['Label']
	ROIList[i]['Mass']   = c_ext_ROI[i]['Mass']
	ROIList[i]['cX']     = c_ext_ROI[i]['cX']
	ROIList[i]['cY']     = c_ext_ROI[i]['cY']
	ROIList[i]['cZ']     = c_ext_ROI[i]['cZ']

    return ROIList[ROIList['Mass']>dust]


def get_blobs(binary_edge_image, mask=1):
    """

    labeled_edge_image, groups = get_blobs(binary_edge_image)

    get the total number of blobs in a 2D or 3D image and convert the 
    binary image (or volume) to labelled regions

    Parameters 
    ----------

    binary_edge_image : {nd_array}
        an binary image/volume

    mask : {int}
        the size of the 2D or 3D connectivity mask. For 2D this is 1, 4 or 8.
	For 3D this is 1, 6, 14 or 28. Mask = 1 is ANY connection in 3x3
	or 3x3x3 mask for 2D or 3D, respectively.

    Returns 
    ----------

    label_image : {nd_array}
        an image/volume with labeled regions from get_blobs() method

    groups : {int}
        number of blobs in image determined by get_blobs() method

    """

    dimensions = binary_edge_image.ndim
    if dimensions == 2:  
        if mask != 1 and mask != 4 and mask != 8:
	    mask = 1 
        [rows, cols] = binary_edge_image.shape
        labeled_edge_image_or_vol = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    elif dimensions == 3:
        if mask != 1 and mask != 6 and mask != 14 and mask != 28:
	    mask = 1 
        [layers, rows, cols] = binary_edge_image.shape
        labeled_edge_image_or_vol = NP.zeros(layers*rows*cols, dtype=NP.uint16).reshape(layers, rows, cols)
    else:
        labeled_edge_image_or_vol = None
	groups = 0
        return labeled_edge_image_or_vol, groups

    groups = S.get_blobs(binary_edge_image, labeled_edge_image_or_vol, mask)

    return labeled_edge_image_or_vol, groups


def sobel_edges(sobel_edge_image, sobel_stats, mode=1, sobel_threshold=0.3):
    """
    sobel_edge = sobel_edges(sobel_edge_image, sobel_stats, mode=1, sobel_threshold=0.3)
    take sobel-filtered image and return binary edges

    Parameters 
    ----------

    sobel_edge_image : {nd_array}
        edge-filtered image from sobel_image() method

    sobel_stats : {dictionary}
        mean and nonzero min, max of sobel filtering

    mode : {0, 1}, optional
        threshold based on histogram mean(0) or mode(1)

    sobel_threshold : {float}, optional
        low threshold applied to edge filtered image for edge generation

    Returns 
    ----------

    sobel_edge : {nd_array}
        binary edge-image

    """
    [rows, cols] = sobel_edge_image.shape
    sobel_edge = NP.zeros(rows*cols, dtype=NP.uint16).reshape(rows, cols)
    S.sobel_edges(sobel_edge_image, sobel_edge, sobel_stats['ave_gt0'], sobel_stats['min_gt0'],
                  sobel_stats['max_gt0'], mode, sobel_threshold)

    return sobel_edge


def sobel_image(filtered_slice):
    """
    sobel_edge_image, sobel_stats = sobel_image(filtered_slice)

    take 2D raw or filtered image slice and get sobel-filtered image 

    Parameters 
    ----------

    filtered_slice : {nd_array}
        raw or pre-processed (filtered and thresholded) 2D image 

    Returns 
    ----------

    sobel_edge_image : {nd_array}
        edge-filtered image from sobel_image() method

    sobel_stats : {dictionary}
        mean and nonzero min, max of sobel filtering

    """
    filtered_slice = filtered_slice.astype(NP.float64)
    [rows, cols] = filtered_slice.shape
    sobel_edge_image = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    pAve, min_value, max_value = S.sobel_image(filtered_slice, sobel_edge_image)

    # can replace this with numpy calls for the image stats. but C ext is faster
    #   S.sobel_image(filtered_slice, sobel_edge_image)
    #   sobel_mask = sobel_edge_image[sobel_edge_image>0]
    #   pAve = sobel_mask.mean()
    #   min_value = sobel_mask.min()
    #   max_value = sobel_mask.max()

    sobel_stats= {'ave_gt0' : pAve, 'min_gt0': min_value, 'max_gt0': max_value} 

    return sobel_edge_image, sobel_stats

def pre_filter(slice, filter, low_threshold=0, high_threshold=0, conv_binary=0):
    """
    edge_filter = pre_filter(slice, filter, low_threshold=0, high_threshold=slice.max(), conv_binary=0)

    take 2D image slice and filter and pre-filter and threshold prior to segmentation

    Parameters 
    ----------
    slice : {nd_array}
        input 2D image. gets cast to int16

    filter : {dictionary}
        2D filter kernel set from build_2d_kernel()

    low_threshold : {int}, optional
        default 0

    high_threshold : {int}, optional
        default max value of source image

    conv_binary : {0, 1}, optional
        flag to convert edge_filter image to binary valued. default 
	is binary conversion off

    Returns 
    ----------
    edge_filter : {nd_array}
        filtered and thresholded image that can be (optional) binary.

    """

    # make sure the input is 16 bits. this is input to edge machine
    # so can handle raw and 8 bit scaled inputs
    if high_threshold==0:
	# default to the maximum value of the image
        high_threshold = slice.max()

    slice = slice.astype(NP.int16)
    [rows, cols] = slice.shape
    edge_image = NP.zeros(rows*cols, dtype=NP.float64).reshape(rows, cols)
    S.edge_prefilter(low_threshold, high_threshold, filter['kernelSize'], filter['kernel'],
		     slice, edge_image)

    if conv_binary == 1:
        edge_image[edge_image>0] = 1
        edge_image = edge_image.astype(NP.uint16)

    return edge_image

def get_max_bounding_box(ROI):
    """
    bounding_box = get_max_bounding_box(ROI)

    take an ROI structure and find the maximum area bounding box

    Parameters 
    ----------

    ROI : {dictionary}
        the ROI is the automatically extracted blob regions of interest
	and contains the rectangular bounding box of each blob.

    Returns 
    ----------

    bounding_box : {dictionary}
        the Left, Right, Top Bottom and Label of the LARGEST bounding box in the ROI

    """
    max_index = ROI[:]['Mass'].argmax()
    bounding_box = {'Left' : ROI[max_index]['Left'], 'Right' : ROI[max_index]['Right'],
		    'Top' : ROI[max_index]['Top'], 'Bottom' : ROI[max_index]['Bottom'],
		    'Label' : ROI[max_index]['Label']} 

    return bounding_box 

def get_all_bounding_boxes(ROI):
    """
    measures = get_all_bounding_boxes(ROI)

    get all bounding boxes in the ROI (feature) dictionary

    Parameters 
    ----------

    ROI : {dictionary}
        the ROI is the automatically extracted blob regions of interest
	and contains the rectangular bounding box of each blob.

    Returns 
    ----------

    measures : {dictionary}
        the Left, Right, Top and Bottom of all bounding boxes in the ROI

    """

    number = ROI.size
    indices = range(0, ROI.size)
    _shortstruct = NP.dtype([('left', 'i'),
                             ('right', 'i'),
                             ('top', 'i'),
                             ('bottom', 'i')])
    measures = NP.zeros(number, dtype=_shortstruct)
    for i in indices:
	measures[i]['left']   = ROI[i]['Left']
	measures[i]['right']  = ROI[i]['Right']
	measures[i]['top']    = ROI[i]['Top']
	measures[i]['bottom'] = ROI[i]['Bottom']

    return measures


def build_2d_kernel(aperature=21, hiFilterCutoff=10.0):
    """

    FIRFilter = build_2d_kernel(aperature, hiFilterCutoff)

    build hamming-windowed FIR filter with sinc kernel

    Parameters 
    ----------

    aperature : {int}, optional
        the number of coefficients in the filter. default is 21. needs to be ODD

    hiFilterCutoff : {float}
        the upper cutoff in digital frequency units

    Returns 
    ----------

    FIRFilter : {dictionary}
        filter kernel

    """

    rad = math.pi / 180.0
    HalfFilterTaps = (aperature-1) / 2
    kernel = NP.zeros((aperature), dtype=NP.float64)
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
    DGFilter = build_d_gauss_kernel(gWidth, sigma)

    build the derivative of Gaussian kernel for Canny edge filter

    Parameters 
    ----------
    gWdith : {int}, optional
         width of derivative of Gaussian kernel.
	 default value is 20

    sigma : {float}, optional
        sigma term of derivative of Gaussian kernel
	 default value is 1.0

    Returns 
    ----------

    DGFilter : {dictionary}
        filter kernel

    """

    kernel  = NP.zeros((1+gWidth), dtype=NP.float64)
    indices = range(0, gWidth)  

    for i in indices:
        kernel[i]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
        kernel[i] *= -(i / (sigma * sigma))

    DGFilter= {'kernelSize' : gWidth, 'coefficients': kernel} 

    return DGFilter

def build_morpho_thin_masks():

    """
    MATFilter = build_morpho_thin_masks()

    build 2 sets (J and K) of 8 3x3 morphology masks (structuring elements)
    to implement thinning (medial axis transformation - MAT)


    Parameters 
    ----------

    None

    Returns 
    ----------

    MATFilter : {dictionary}
        morphology filter kernels. there are 2 sets of 8 3x3 masks

    """

    # (layers, rows, cols)
    size = (8*3*3)
    J_mask = NP.zeros(size, dtype=NP.int16)
    K_mask = NP.zeros(size, dtype=NP.int16)

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


def build_laws_kernel():

    """
    LAWSFilter = build_laws_kernel()

    build 6 length-7 Law's texture filter masks
    mask names are: 'L', 'S', 'E', 'W', 'R', 'O'


    Parameters 
    ----------

    None

    Returns 
    ----------

    LAWSFilter : {dictionary}
        a set of 6 length-7 Laws texture kernels

    """
    aperature = (6, 7)
    coefficients = NP.zeros((aperature), dtype=NP.float64)
    names = ('L', 'E', 'S', 'W', 'R', 'O' )

    coefficients[0, :] =  ( 1.0,  6.0,  15.0, 20.0,  15.0,  6.0,  1.0)
    coefficients[1, :] =  (-1.0, -4.0,  -5.0,  0.0,   5.0,  4.0,  1.0)
    coefficients[2, :] =  (-1.0, -2.0,   1.0,  4.0,   1.0, -2.0, -1.0)
    coefficients[3, :] =  (-1.0,  0.0,   3.0,  0.0,  -3.0,  0.0,  1.0)
    coefficients[4, :] =  ( 1.0, -2.0,  -1.0,  4.0,  -1.0, -2.0,  1.0)
    coefficients[5, :] =  (-1.0,  6.0, -15.0, 20.0, -15.0,  6.0, -1.0)

    LAWSFilter= {'numKernels' : 6, 'kernelSize' : 7, 'filters' : 21,
		 'coefficients': coefficients, 'names': names} 

    return LAWSFilter

def build_laws_masks(LAWSFilter):

    """
    masks = build_laws_masks(LAWSFilter)

    takes the Laws Filter dictionary and builds the 21 7x7 kernel masks that are
    used in Laws texture feature extraction. 

    Parameters 
    ----------

    LAWSFilter : {dictionary}
        a set of 6 length-7 Laws texture kernels

    Returns 
    ----------

    masks : {list}
        a list of 21 7x7 kernels (2D nd_array)

    Examples:
    use this along with FFT Pack to view the spatial frequency response. Create a 256x256 zero 
    array and pad the first 7x7 with the Laws kernel and then get the 2D power spectrum and display
    with pylab


    LAWSFilter = build_laws_kernel()
    mask = build_laws_masks(LAWSFilter)

    mask_2 = masks[2]

    z = NP.zeros(256*256, dtype=NP.float32).reshape(256, 256)
    z[0:7, 0:7] = mask_0
    x = abs(fftshift(fft2(z)))
    pylab.imshow(x)

    """

    outer_indices = range(0, LAWSFilter['numKernels'])
    mask_array = {}
    count = 0
    for i in outer_indices:
	rowFilter = LAWSFilter['coefficients'][i]
	colFilter = LAWSFilter['coefficients'][i]
	matrix = NP.outer(rowFilter, colFilter)
	mask_array[count] = 2.0*matrix
	count = count + 1 
        inner_indices = range(i+1, LAWSFilter['numKernels'])
        for j in inner_indices:
	    colFilter = LAWSFilter['coefficients'][j]
	    matrix = NP.outer(rowFilter, colFilter) + NP.outer(colFilter, rowFilter)
	    mask_array[count] = matrix
	    count = count + 1 

    return mask_array


#
#    test pattern generators for demo and test
#

def build_test_texture_discs():
    """
    discs = build_test_texture_discs()
    
    builds 4 discs with plane wave texture. used for test and demo

    Parameters 
    ----------

    None

    Returns 
    ----------

    discs : {nd_array}
        a 512x512 image with 4 test discs (one per quadrant)

    """
    rows = 512
    cols = 512
    rad  = NP.pi / 180.0
    test_image = NP.zeros(rows*cols, dtype=NP.float32).reshape(rows, cols)
    [a, b] = NP.mgrid[0:rows, 0:cols]

    test_image[0:255,0:255]     = NP.sin(4.0*rad*a[0:255,0:255])  + NP.sin(-4.0*rad*b[0:255,0:255]) 
    test_image[256:511,256:511] = NP.sin(24.0*rad*a[0:255,0:255]) + NP.sin(20.0*rad*b[0:255,0:255])

    test_image = test_image + test_image.min()
    discs = build_test_unit_discs()
    discs = discs * test_image 

    return discs


def build_test_discs():
    """
    test_image = build_test_discs()
    build 4 discs of equal radius and different mean values for edge/blob testing
    
    Parameters 
    ----------

    None

    Returns 
    ----------

    test_image : {nd_array}
        a 512x512 image with 4 test discs (one per quadrant)

    """
    radius = 50
    rows   = 512
    cols   = 512
    test_image = NP.zeros(rows*cols, dtype=NP.int16).reshape(rows, cols)
    y_indices = NP.array(range(-radius, radius+1))
    center_x = rows / 4
    center_y = cols / 4

    for i in y_indices:
	x = math.sqrt(float(radius)**2 - float(i)**2)
	# different raw mean levels
	test_image[1*center_y+i, 1*center_x-x:1*center_x+x] = 80
	test_image[1*center_y+i, 3*center_x-x:3*center_x+x] = 90
	test_image[3*center_y+i, 1*center_x-x:1*center_x+x] = 100
	test_image[3*center_y+i, 3*center_x-x:3*center_x+x] = 110

    return test_image


def build_test_unit_discs():
    """
    test_image = build_test_unit_discs()
    build 2 discs of equal radius and same mean values for texture testing
    
    Parameters 
    ----------

    None

    Returns 
    ----------

    test_image : {nd_array}
        a 512x512 image with 4 test discs (one per quadrant)

    """
    radius = 50
    rows   = 512
    cols   = 512
    test_image = NP.zeros(rows*cols, dtype=NP.int16).reshape(rows, cols)
    y_indices = NP.array(range(-radius, radius+1))
    center_x = rows / 4
    center_y = cols / 4

    for i in y_indices:
	x = math.sqrt(float(radius)**2 - float(i)**2)
	# different raw mean levels
	test_image[1*center_y+i, 1*center_x-x:1*center_x+x] = 100
	test_image[3*center_y+i, 3*center_x-x:3*center_x+x] = 100

    return test_image


def build_test_impulses():
    """
    test_image = build_test_impulses()

    build 4 test impulses discs centered in the 4 discs. used
    for testing filter kernels, esp. Laws' filter kernel. Filtering
    with these test patterns will return Law's kernel outer product matrices.

    Parameters 
    ----------

    None

    Returns 
    ----------

    test_image : {nd_array}
        a 512x512 image with 4 test discs (one per quadrant)

    """
    rows = 512
    cols = 512
    test_image = NP.zeros(rows*cols, dtype=NP.int16).reshape(rows, cols)
    test_image[128,128] = 1 
    test_image[378,128] = 1 
    test_image[128,378] = 1 
    test_image[378,378] = 1 

    return test_image



