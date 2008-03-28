import math
import os
import numpy as NP
import scipy.ndimage._register as R
import scipy.special  as SP
import scipy.ndimage  as NDI
import scipy.optimize as OPT
import time
import glob

# Issue warning regarding heavy development status of this module
import warnings
_msg = "The registration code is under heavy development and therefore the \
public API will change in the future.  The NIPY group is actively working on \
this code, and has every intention of generalizing this for the Scipy \
community.  Use this module minimally, if at all, until it this warning is \
removed."
warnings.warn(_msg, UserWarning)

# TODO:  Add docstrings for public functions in extension code.
# Add docstrings to extension code.
#from numpy.lib import add_newdoc
#add_newdoc('scipy.ndimage._register', 'register_histogram',
#    """A joint histogram used for registration module.
#    """)


#
#  ---- co-registration and IO  ---- 
#

def resize_image(imageG, imageF_mat):
	"""
	zoom_image = resize_image(source_image, reference_image[mat])

	Fractional resample source_image to reference_imagesize. The
	resample is implemented with 3D cubic spline. The reference
	image [mat] is the 4x4 voxel-to-physical conversion matrix.
	
	Parameters 
	----------

	imageG : {dictionary} 
	    imageG is the source image to be resized. it is a dictionary with
	    the data as an ndarray in the ['data'] component.

	reference_image[mat] : {ndarray}
	    refernce_image is the image whose sampling dimensions the source
	    image is to be remapped to. [mat] refers to the component
	    of the image dictionary, reference_image['mat'] that is the
	    sampling dimensions.

	Returns 
	-------
	zoom_image : {dictionary}

	Examples
	--------

	>>> import _registration as reg
	>>> measures, imageF_anat, fmri_series = reg.demo_MRI_coregistration()

	>>> resampled_fmri = reg.resize_image(fmri_series[10], imageF_anat['mat'])

	image 10 in the fmri_series is resampled to imageF_anat coordinates

	"""

	Z = NP.zeros(3, dtype=NP.float64);
	# get the zoom
	Z[0] = imageG['mat'][0][0] / imageF_mat[0][0]
	Z[1] = imageG['mat'][1][1] / imageF_mat[1][1]
	Z[2] = imageG['mat'][2][2] / imageF_mat[2][2]

	# new volume dimensions (rounded)
	D = NP.zeros(3, dtype=NP.int32);
	D[0] = int(float(imageG['dim'][0])*Z[0]+0.5)
	D[1] = int(float(imageG['dim'][1])*Z[1]+0.5)
	D[2] = int(float(imageG['dim'][2])*Z[2]+0.5)

	M = NP.eye(4, dtype=NP.float64);
	# for the test data, set the xyz voxel sizes for fMRI volume
	M[0][0] = imageG['mat'][0][0]/Z[0]
	M[1][1] = imageG['mat'][1][1]/Z[1]
	M[2][2] = imageG['mat'][2][2]/Z[2]

    	image = NP.zeros(D[2]*D[1]*D[0], dtype=NP.uint8).reshape(D[2], D[0], D[1])
	mode  = 2
	scale = 0
	R.register_volume_resample(imageG['data'], image, Z, scale, mode)
	F = NP.zeros(3, dtype=NP.float64);
	zoom_image = {'data' : image, 'mat' : M, 'dim' : D, 'fwhm' : F}

	return zoom_image

def remap_image(image, parm_vector, resample='linear'):
	"""
	remaped_image = remap_image(image, parm_vector, resample='linear')

	rotates and translates image using the 3 angles and 3 translations in the 6-dim
	parm_vector. The mapping is stored in the 4x4 M_inverse matrix from the get_inverse_mapping
	method.

	Parameters 
	----------
	image : {dictionary} 
	    image is the source image to be remapped. it is a dictionary with
	    the data as an ndarray in the ['data'] component.

	parm_vector : {ndarray}
	    parm_vector is the 6-dimensional vector (3 angles, 3 translations)
	    generated from the registration.

	resample : {'linear', 'cubic'}, optional


	Returns 
	-------
	remaped_image : {dictionary}

	Examples
	--------
	    image = fmri_series[i]
	    x[0:6] = measures[i]['align_rotate'][0:6]
	    # overwrite the fMRI volume with the aligned volume
	    fmri_series[i] = remap_image(image, x, resample='cubic')

	"""

	#
	# remap imageG to coordinates of imageF (creates imageG')
	# use the 6 dim parm_vector (3 angles, 3 translations) to remap
	#
	M_inverse = get_inverse_mappings(parm_vector)
	(layers, rows, cols) = image['data'].shape
	# allocate the zero image
	remaped_image = NP.zeros(layers*rows*cols, dtype=NP.uint8).reshape(layers, rows, cols)
	remaped_image = {'data' : remaped_image, 'mat' : image['mat'], 
			 'dim' : image['dim'], 'fwhm' : image['fwhm']}
	imdata = build_structs()

	if resample == 'linear':
	    # trilinear interpolation mapping.
	    R.register_linear_resample(image['data'], remaped_image['data'], M_inverse, imdata['step'])
	elif resample == 'cubic':
	    # tricubic convolve interpolation mapping. 
	    R.register_cubic_resample(image['data'], remaped_image['data'], M_inverse, imdata['step'])

	return remaped_image

def get_inverse_mappings(parm_vector):
	"""
	M_inverse = get_inverse_mappings(parm_vector)

	takes the 6-dimensional rotation/translation vector and builds the inverse
	4x4 mapping matrix M_inverse that will map imageG to imageF orientation

	Parameters 
	----------
	parm_vector : {nd_array} 

	Returns 
	-------
	M_inverse : {nd_array}

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> array = NP.zeros(6, dtype=float)
	>>> M = reg.get_inverse_mappings(array)
	>>> M 

	array([
	[ 1.,  0.,  0.,  0.],
       	[ 0.,  1.,  0.,  0.],
       	[ 0.,  0.,  1.,  0.],
       	[ 0.,  0.,  0.,  1.]])

	"""
	# get the inverse mapping to rotate the G matrix to F space following registration
	imdata = build_structs()
	# inverse angles and translations
	imdata['parms'][0] = -parm_vector[0]
	imdata['parms'][1] = -parm_vector[1]
	imdata['parms'][2] = -parm_vector[2]
	imdata['parms'][3] = -parm_vector[3]
	imdata['parms'][4] = -parm_vector[4]
	imdata['parms'][5] = -parm_vector[5]
	M_inverse = build_rotate_matrix(imdata['parms'])
	return M_inverse

def python_coreg(image1, image2, imdata, ftype=1, smimage=0, lite=0, smhist=0,
		 method='nmi', opt_method='powell'):
	"""
	parm_vector = python_coreg(image1, image2, imdata, ftype=1, smimage=0, lite=0,
				   smhist=0, method='nmi', opt_method='powell'):

	takes two images and the image data descriptor (imdata) and determines the optimal 
	alignment of the two images (measured by mutual information or cross correlation) 
	using optimization search of 3 angle and 3 translation parameters. The optimization 
	uses either the Powell or Conjugate Gradient methods in the scipy optimization 
	package. The optimal parameter is returned.

	Parameters 
	----------
	image1 : {dictionary} 
	    image1 is the source image to be remapped during the registration. 
	    it is a dictionary with the data as an ndarray in the ['data'] component.
	image2 : {dictionary} 
	    image2 is the reference image that image1 gets mapped to. 
	imdata : {dictionary} 
	    image sampling and optimization information.
	ftype : {0, 1}, optional
	    flag for type of low pass filter. 0 is Gauss-Spline
	    1 is pure Gauss. Sigma determined from volume sampling info.
	smimage : {0, 1}, optional
	    flag for volume 3D low pass filtering of image 2.
	    0 for no filter, 1 for do filter.
	lite : {0, 1}, optional
	    lite of 1 is to jitter both images during resampling. 0
	    is to not jitter. jittering is for non-aliased volumes.
	smhist: {0, 1}, optional
	    flag for joint histogram low pass filtering. 0 for no filter,
	    1 for do filter.
	method: {'nmi', 'mi', 'ncc', 'ecc'}, optional
	    flag for type of registration metric. nmi is normalized mutual
	    information; mi is mutual information; ecc is entropy cross
	    correlation; ncc is normalized cross correlation.
	opt_method: {'powell', 'hybrid'}, optional
	    registration is two pass. Pass 1 is low res to get close to alignment
	    and pass 2 starts at the pass 1 optimal alignment. In powell pass 1 and
	    2 are powell, in hybrid pass 2 is conjugate gradient.


	Returns 
	-------
	parm_vector : {nd_array}
	    this is the optimal alignment (6-dim) array with 3 angles and
	    3 translations.

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg

	>>> image1, image2, imdata = reg.demo_MRI_volume_align()
	>>> parm_vector = python_coreg(image1, image2, imdata)

	"""
    	start = time.time()
	# smooth of the images
	if smimage: 
	    image_F_xyz2 = filter_image_3D(image2['data'], image2['fwhm'], ftype)
	    image2['data'] = image_F_xyz2
	parm_vector = multires_registration(image1, image2, imdata, lite, smhist, method, opt_method)
    	stop = time.time()
	print 'Total Optimizer Time is ', (stop-start)
	return parm_vector

def multires_registration(image1, image2, imdata, lite, smhist, method, opt_method):
	"""
	x = multires_registration(image1, image2, imdata, lite, smhist, method, opt_method)

	to be called by python_coreg() which optionally does 3D image filtering and 
	provies timing for registration.

	Parameters 
	----------

	image1 : {dictionary} 
	    image1 is the source image to be remapped during the registration. 
	    it is a dictionary with the data as an ndarray in the ['data'] component.
	image2 : {dictionary} 
	    image2 is the reference image that image1 gets mapped to. 
	imdata : {dictionary} 
	    image sampling and optimization information.
	lite : {integer}
	    lite of 1 is to jitter both images during resampling. 0
	    is to not jitter. jittering is for non-aliased volumes.
	smhist: {integer}
	    flag for joint histogram low pass filtering. 0 for no filter,
	    1 for do filter.
	method: {'nmi', 'mi', 'ncc', 'ecc'}
	    flag for type of registration metric. nmi is normalized mutual
	    information; mi is mutual information; ecc is entropy cross
	    correlation; ncc is normalized cross correlation.
	opt_method: {'powell', 'hybrid'}
	    registration is two pass. Pass 1 is low res to get close to alignment
	    and pass 2 starts at the pass 1 optimal alignment. In powell pass 1 and
	    2 are powell, in hybrid pass 2 is conjugate gradient.

	Returns 
	-------
	x : {nd_array}
	    this is the optimal alignment (6-dim) array with 3 angles and
	    3 translations.

	Examples
	--------

	(calling this from python_coreg which optionally filters image2)
	>>> import numpy as NP
	>>> import _registration as reg
	>>> image1, image2, imdata = reg.demo_MRI_volume_align()
	>>> parm_vector = python_coreg(image1, image2, imdata)

	"""
	ret_histo=0
	# zero out the start parameter; but this may be set to large values 
	# if the head is out of range and well off the optimal alignment skirt
	imdata['parms'][0:5] = 0.0
	# make the step a scalar to can put in a multi-res loop
	loop = range(imdata['sample'].size)
    	x = imdata['parms']
	for i in loop:
	    step = imdata['sample'][i]
	    imdata['step'][:] = step
	    optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist,
			    method, ret_histo)
	    p_args = (optfunc_args,)
	    if opt_method=='powell':
		print 'POWELL multi-res registration step size ', step
		print 'vector ', x
    	        x = OPT.fmin_powell(optimize_function, x, args=p_args,
				    callback=callback_powell) 
	    elif opt_method=='cg':
		print 'CG multi-res registration step size ', step
		print 'vector ', x
    	        x = OPT.fmin_cg(optimize_function, x, args=p_args, callback=callback_cg) 
	    elif opt_method=='hybrid':
		if i==0:
		    print 'Hybrid POWELL multi-res registration step size ', step
		    print 'vector ', x
		    lite = 0
	    	    optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist,
				    method, ret_histo)
	    	    p_args = (optfunc_args,)
    	            x = OPT.fmin_powell(optimize_function, x, args=p_args, callback=callback_powell) 
	        elif i==1:
		    print 'Hybrid CG multi-res registration step size ', step
		    print 'vector ', x
		    lite = 1
	    	    optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, 
				    smhist, method, ret_histo)
	    	    p_args = (optfunc_args,)
    	            x = OPT.fmin_cg(optimize_function, x, args=p_args, callback=callback_cg) 

	return x


def callback_powell(x):
	"""
	called by optimize.powell only. prints current parameter vector.
	"""
	print 'Parameter Vector from Powell: - '
	print x
	return

def callback_cg(x):
	"""
	called by optimize.cg only. prints current parameter vector.
	"""
	print 'Parameter Vector from Conjugate Gradient: - '
	print x
	return

def smooth_kernel(fwhm, x, ktype=1):
	"""
	kernel = smooth_kernel(fwhm, x, ktype=1)

	smooth_kernel creates filter kernel based on image sampling parameters.
	provide domain of kernel and sampling parameters. 

	Parameters 
	----------
	fwhm : {int}
	    used for kernel width
	x : {nd_array}
	    domain of kernel
	ktype: {1, 2}, optional
	    kernel type. 1 is Gauss convoled with spline, 2 is Gauss


	Returns 
	-------
	kernel : {nd_array}

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> fwhm = 3
	>>> ftype = 2
	>>> p = NP.ceil(2*fwhm).astype(int)
	>>> x = NP.array(range(-p, p+1))
	>>> kernel = reg.smooth_kernel(fwhm, x, ktype=ftype)
	>>> kernel

	array([
	 4.77853772e-06,   1.41575516e-04,   2.26516955e-03,
         1.95718875e-02,   9.13238336e-02,   2.30120330e-01,
         3.13144850e-01,   2.30120330e-01,   9.13238336e-02,
         1.95718875e-02,   2.26516955e-03,   1.41575516e-04,
         4.77853772e-06])

	"""
	eps = 0.00001
	s   = NP.square((fwhm/math.sqrt(8.0*math.log(2.0)))) + eps
	if ktype==1:
	    # from SPM: Gauss kernel convolved with 1st degree B spline
	    w1 = 0.5 * math.sqrt(2.0/s)
	    w2 = -0.5 / s
	    w3 = math.sqrt((s*math.pi) /2.0)
	    kernel = 0.5*(SP.erf(w1*(x+1))*(x+1)       + SP.erf(w1*(x-1))*(x-1)    - 2.0*SP.erf(w1*(x))*(x) + 
	 	          w3*(NP.exp(w2*NP.square(x+1))) + NP.exp(w2*(NP.square(x-1))) - 2.0*NP.exp(w2*NP.square(x)))
	    kernel[kernel<0] = 0
	    kernel = kernel / kernel.sum()  
	else:
	    # Gauss kernel 
	    kernel = (1.0/math.sqrt(2.0*math.pi*s)) * NP.exp(-NP.square(x)/(2.0*s)) 
	    kernel = kernel / kernel.sum()  

	return kernel

def filter_image_3D(imageRaw, fwhm, ftype=2):
	"""
	image_F_xyz = filter_image_3D(imageRaw, fwhm, ftype=2):
	does 3D separable digital filtering using scipy.ndimage.correlate1d

	Parameters 
	----------
	imageRaw : {nd_array}
	    the unfiltered 3D volume image
	fwhm : {int}
	    used for kernel width
	ktype: {1, 2}, optional
	    kernel type. 1 is Gauss convoled with spline, 2 is Gauss

	Returns 
	-------
	image_F_xyz : {nd_array}
	    3D filtered volume image

	Examples
	--------

	>>> import _registration as reg
	>>> image1, image2, imdata = reg.demo_MRI_volume_align()
	>>> ftype = 1
	>>> image_Filter_xyz = filter_image_3D(image1['data'], image1['fwhm'], ftype)
	>>> image1['data'] = image_Filter_xyz
	"""

	p = NP.ceil(2*fwhm[0]).astype(int)
	x = NP.array(range(-p, p+1))
	kernel_x = smooth_kernel(fwhm[0], x, ktype=ftype)
	p = NP.ceil(2*fwhm[1]).astype(int)
	x = NP.array(range(-p, p+1))
	kernel_y = smooth_kernel(fwhm[1], x, ktype=ftype)
	p = NP.ceil(2*fwhm[2]).astype(int)
	x = NP.array(range(-p, p+1))
	kernel_z = smooth_kernel(fwhm[2], x, ktype=ftype)
	output=None
	# 3D filter in 3 1D separable stages
	axis = 0
	image_F_x   = NDI.correlate1d(imageRaw,   kernel_x, axis, output)
	axis = 1
	image_F_xy  = NDI.correlate1d(image_F_x,  kernel_y, axis, output)
	axis = 2
	image_F_xyz = NDI.correlate1d(image_F_xy, kernel_z, axis, output)
	return image_F_xyz  

def build_fwhm(M, S):
	"""
	fwhm = build_fwhm(M, S)

	builds the low pass filter kernel sigma value from the image pixel sampling

	Parameters 
	----------
	M : {nd_array}
	    input 4x4 voxel to physical map matrix (called 'MAT')

	S : {nd_array}
	    1x3 sample increment array. should be = (1, 1, 1)

	Returns 
	-------
	fwhm : {nd_array}
	    the 3D Gaussian kernel width

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> anat_desc = reg.load_anatMRI_desc()
	>>> image1 = reg.load_volume(anat_desc, imagename='ANAT1_V0001.img')
	>>> imdata = reg.build_structs()
	>>> image1['fwhm'] = reg.build_fwhm(image1['mat'], imdata['step'])

	"""
	view_3x3 = NP.square(M[0:3, 0:3])
	# sum the elements inn the first row
	vxg = NP.sqrt(view_3x3.sum(axis=0))
	# assumes that sampling is the same for xyz
 	size = NP.array([1,1,1])*S[0]
	x = NP.square(size) - NP.square(vxg)
	# clip
	x[x<0] = 0
	fwhm = NP.sqrt(x) / vxg
	# pathology when stepsize = 1 for MAT equal to the identity matrix
	fwhm[fwhm==0] = 1
	# return the 3D Gaussian kernel width (xyz)
	return fwhm 

def optimize_function(x, optfunc_args):
	"""
	cost = optimize_function(x, optfunc_args)    --- OR ---
	cost, joint_histogram = optimize_function(x, optfunc_args)   

 	computes the alignment between 2 volumes using cross correlation or mutual
	information metrics. In both the 8 bit joint histogram of the 2 images is
	computed. The 8 bit scaling is done using an integrated histogram method and
	is called prior to this.

	Parameters 
	----------
	x : {nd_array}
	    this is the current (6-dim) array with 3 angles and 3 translations.

	optfunc_args : {tuple}
	    this is a tuple of 8 elements that is formed in the scipy.optimize powell
	    and cg (conjugate gradient) functions. this is why the elements need to be
	    a tuple. The elements of optfunc_args are: 

	    image_F       : {dictionary} 
	        image_F is the source image to be remapped during the registration. 
	        it is a dictionary with the data as an ndarray in the ['data'] component.
	    image_G       : {dictionary} 
	        image_G is the reference image that image_F gets mapped to. 
	    sample_vector : {nd_array} 
	        sample in x,y,x. should be (1,1,1)
	    fwhm          : {nd_array} 
	        Gaussian sigma
	    do_lite       : {0, 1} 
	        lite of 1 is to jitter both images during resampling. 
	        0 is to not jitter. jittering is for non-aliased volumes.
	    smooth        : {0, 1} 
	        flag for joint histogram low pass filtering. 0 for no filter,
	        1 for do filter.
	    method        : {'nmi', 'mi', 'ncc', 'ecc', 'mse'}
	        flag for type of registration metric. nmi is normalized mutual
	        information; mi is mutual information; ecc is entropy cross
	        correlation; ncc is normalized cross correlation. mse is mean
		square error. with mse there is no joint histogram.
	    ret_histo     : {0, 1} 
	        if 0 return is: cost 
	        if 0 return is: cost, joint_histogram  

	Returns 
	-------
	    cost : {float}
	        the negative of one of the mutual information metrics
		or negative cross correlation. use negative as the optimization
		is minimization.

	    --- OR --- (if ret_histo = 1)

	    cost : {float}
	        the negative of one of the mutual information metrics
		or negative cross correlation. use negative as the optimization
		is minimization.

	    joint_histogram : {nd_array}
	        the 2D (256x256) joint histogram of the two volumes


	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> anat_desc = reg.load_anatMRI_desc()
	>>> image1 = reg.load_volume(anat_desc, imagename='ANAT1_V0001.img')
	>>> image2 = reg.load_volume(anat_desc, imagename='ANAT1_V0001.img')
	>>> imdata = reg.build_structs()
	>>> image1['fwhm'] = reg.build_fwhm(image1['mat'], imdata['step'])
	>>> image2['fwhm'] = reg.build_fwhm(image2['mat'], imdata['step'])
	>>> method = 'ncc'
	>>> lite = 1
	>>> smhist = 0
	>>> ret_histo = 1
	>>> optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)
	>>> x = NP.zeros(6, dtype=NP.float64)
	>>> return cost, joint_histogram = reg.optimize_function(x, optfunc_args)


	"""

	image_F       = optfunc_args[0]
	image_G       = optfunc_args[1]
	sample_vector = optfunc_args[2]
	fwhm          = optfunc_args[3]
	do_lite       = optfunc_args[4]
	smooth        = optfunc_args[5]
	method        = optfunc_args[6]
	ret_histo     = optfunc_args[7]

	rot_matrix = build_rotate_matrix(x)
	cost = 0.0
	epsilon = 2.2e-16 
	# image_G is base image
	# image_F is the to-be-rotated image
	# rot_matrix is the 4x4 constructed (current angles and translates) transform matrix
	# sample_vector is the subsample vector for x-y-z

	F_inv = NP.linalg.inv(image_F['mat'])
	composite = NP.dot(F_inv, image_G['mat'])
	composite = NP.dot(composite, rot_matrix)

	if method == 'mse':
	    #
	    # mean squard error method
	    #

	    (layers, rows, cols) = image_F['data'].shape
	    # allocate the zero image
	    remap_image_F = NP.zeros(layers*rows*cols, dtype=NP.uint8).reshape(layers, rows, cols)
	    imdata = build_structs()
	    # trilinear interpolation mapping.
	    R.register_linear_resample(image_F['data'], remap_image_F, composite,
			               imdata['step'])
	    cost = (NP.square(image_G['data']-remap_image_F)).mean()

	    return cost

	else:
	    #
	    # histogram-based methods (nmi, ncc, mi, ecc)
	    #

	    # allocate memory for 2D histogram
	    joint_histogram = NP.zeros([256, 256], dtype=NP.float64);

	    if do_lite: 
	        R.register_histogram_lite(image_F['data'], image_G['data'], composite,
			                  sample_vector, joint_histogram)
	    else:
	        R.register_histogram(image_F['data'], image_G['data'], composite,
			             sample_vector, joint_histogram)

	    # smooth the histogram
	    if smooth: 
	        p = NP.ceil(2*fwhm[0]).astype(int)
	        x = NP.array(range(-p, p+1))
	        kernel1 = smooth_kernel(fwhm[0], x)
	        p = NP.ceil(2*fwhm[1]).astype(int)
	        x = NP.array(range(-p, p+1))
	        kernel2 = smooth_kernel(fwhm[1], x)
	        output=None
	        # 2D filter in 1D separable stages
	        axis = 0
	        result = NDI.correlate1d(joint_histogram, kernel1, axis, output)
	        axis = 1
	        joint_histogram = NDI.correlate1d(result, kernel1, axis, output)

	    joint_histogram += epsilon # prevent log(0) 
	    # normalize the joint histogram
	    joint_histogram /= joint_histogram.sum() 
	    # get the marginals
	    marginal_col = joint_histogram.sum(axis=0)
	    marginal_row = joint_histogram.sum(axis=1)

	    if method == 'mi':
	        # mutual information
	        marginal_outer = NP.outer(marginal_col, marginal_row)
	        H = joint_histogram * NP.log(joint_histogram / marginal_outer)  
	        mutual_information = H.sum()
	        cost = -mutual_information

	    elif method == 'ecc':
	        # entropy correlation coefficient 
	        marginal_outer = NP.outer(marginal_col, marginal_row)
	        H = joint_histogram * NP.log(joint_histogram / marginal_outer)  
	        mutual_information = H.sum()
	        row_entropy = marginal_row * NP.log(marginal_row)
	        col_entropy = marginal_col * NP.log(marginal_col)
	        ecc  = -2.0*mutual_information/(row_entropy.sum() + col_entropy.sum())
	        cost = -ecc

	    elif method == 'nmi':
	        # normalized mutual information
	        row_entropy = marginal_row * NP.log(marginal_row)
	        col_entropy = marginal_col * NP.log(marginal_col)
	        H = joint_histogram * NP.log(joint_histogram)  
	        nmi = (row_entropy.sum() + col_entropy.sum()) / (H.sum())
	        cost = -nmi

	    elif method == 'ncc':
	        # cross correlation from the joint histogram 
	        r, c = joint_histogram.shape
	        i = NP.array(range(1,c+1))
	        j = NP.array(range(1,r+1))
	        m1 = (marginal_row * i).sum()
	        m2 = (marginal_col * j).sum()
	        sig1 = NP.sqrt((marginal_row*(NP.square(i-m1))).sum())
	        sig2 = NP.sqrt((marginal_col*(NP.square(j-m2))).sum())
	        [a, b] = NP.mgrid[1:c+1, 1:r+1]
	        a = a - m1
	        b = b - m2
	        # element multiplies in the joint histogram and grids
	        H = ((joint_histogram * a) * b).sum()
	        ncc = H / (NP.dot(sig1, sig2)) 
	        cost = -ncc

	    if ret_histo:
	        return cost, joint_histogram 
    	    else:
	        return cost


def build_structs(step=1):
	"""
	img_data = build_structs(step=1)

	builds the image data (imdata) dictionary for later use as parameter
	storage in the co-registration.

	Parameters 
	----------
	step : {int} : optional
	default is 1 and is the sample increment in voxels. This sets the sample
	for x,y,z and is the same value in all 3 axes. only change the default for debug.

	Returns 
	-------
	img_data : {dictionary}

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> imdata = reg.build_structs()

	"""

	# build image data structures here
	P = NP.zeros(6, dtype=NP.float64);
	T = NP.zeros(6, dtype=NP.float64);
	F = NP.zeros(2, dtype=NP.int32);
	S = NP.ones(3,  dtype=NP.int32);
	sample = NP.zeros(2, dtype=NP.int32);
	S[0] = step
	S[1] = step
	S[2] = step
	# image/histogram smoothing
	F[0] = 3
	F[1] = 3
	# subsample for multiresolution registration
	sample[0] = 4
	sample[1] = 2
	# tolerances for angle (0-2) and translation (3-5)
	T[0] = 0.02 
	T[1] = 0.02 
	T[2] = 0.02 
	T[3] = 0.001 
	T[4] = 0.001 
	T[5] = 0.001 
	# P[0] = alpha <=> pitch. + alpha is moving back in the sagittal plane
	# P[1] = beta  <=> roll.  + beta  is moving right in the coronal plane
	# P[2] = gamma <=> yaw.   + gamma is right turn in the transverse plane
	# P[3] = Tx
	# P[4] = Ty
	# P[5] = Tz
	img_data = {'parms' : P, 'step' : S, 'fwhm' : F, 'tol' : T, 'sample' : sample}
	return img_data


def build_rotate_matrix(img_data_parms):
	"""
	rot_matrix = reg.build_rotate_matrix(img_data_parms)

	takes the 6 element vector (3 angles, 3 translations) and build the 4x4 mapping matrix

	Parameters 
	----------
	img_data_parms : {nd_array}
	    this is the current (6-dim) array with 3 angles and 3 translations.

	Returns 
	-------
	rot_matrix: {nd_array}
	    the 4x4 mapping matrix

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> imdata = reg.build_structs()
	>>> x = NP.zeros(6, dtype=NP.float64)
	>>> M = reg.build_rotate_matrix(x)
	>>> M 
	array([[ 1.,  0.,  0.,  0.],
       	       [ 0.,  1.,  0.,  0.],
       	       [ 0.,  0.,  1.,  0.],
       	       [ 0.,  0.,  0.,  1.]])


	"""

	R1 = NP.zeros([4,4], dtype=NP.float64);
	R2 = NP.zeros([4,4], dtype=NP.float64);
	R3 = NP.zeros([4,4], dtype=NP.float64);
	T  = NP.eye(4, dtype=NP.float64);

	alpha = math.radians(img_data_parms[0])
	beta  = math.radians(img_data_parms[1])
	gamma = math.radians(img_data_parms[2])

	R1[0][0] = 1.0
	R1[1][1] = math.cos(alpha)
	R1[1][2] = math.sin(alpha)
	R1[2][1] = -math.sin(alpha)
	R1[2][2] = math.cos(alpha)
	R1[3][3] = 1.0

	R2[0][0] = math.cos(beta)
	R2[0][2] = math.sin(beta)
	R2[1][1] = 1.0
	R2[2][0] = -math.sin(beta)
	R2[2][2] = math.cos(beta)
	R2[3][3] = 1.0

	R3[0][0] = math.cos(gamma)
	R3[0][1] = math.sin(gamma)
	R3[1][0] = -math.sin(gamma)
	R3[1][1] = math.cos(gamma)
	R3[2][2] = 1.0
	R3[3][3] = 1.0

	T[0][0] = 1.0
	T[1][1] = 1.0
	T[2][2] = 1.0
	T[3][3] = 1.0
	T[0][3] = img_data_parms[3]
	T[1][3] = img_data_parms[4]
	T[2][3] = img_data_parms[5]

	rot_matrix = NP.dot(T, R1);
	rot_matrix = NP.dot(rot_matrix, R2);
	rot_matrix = NP.dot(rot_matrix, R3);

	return rot_matrix


def load_volume(imagedesc, imagename=None, threshold=0.999, debug=0):

	"""
	image = load_volume(imagedesc, imagename=None, threshold=0.999, debug=0)  --- OR ---
	image, h, ih, index = load_volume(imagedesc, imagename=None, threshold=0.999, debug=0)

	gets an image descriptor and optional filename and returns a scaled 8 bit volume. The
	scaling is designed to make full use of the 8 bits (ignoring high amplitude outliers).
	The current method uses numpy fromfile and will be replaced by neuroimage nifti load.

	Parameters 
	----------
	imagedesc : {dictionary} 
	    imagedesc is the descriptor of the image to be read. 

	imagename : {string} : optional
	    name of image file. No name creates a blank image that is used for creating
	    a rotated test image or image rescaling.

	threshold : {float} : optional
	    this is the threshold for upper cutoff in the 8 bit scaling. The volume histogram
	    and integrated histogram is computed and the upper amplitude cutoff is where the 
	    integrated histogram crosses the value set in the threshold. setting threshold to
	    1.0 means the scaling is done over the min to max amplitude range.

	debug : {0, 1} : optional
	    when debug=1 the method returns the volume histogram, integrated histogram and the 
	    amplitude index where the provided threshold occured.

	Returns 
	-------
	image : {dictionary}
	    the volume data assoicated with the filename or a blank volume of the same
	    dimensions as specified in imagedesc.

	--- OR --- (if debug = 1)

	image : {dictionary}
	    the volume data assoicated with the filename or a blank volume of the same
	    dimensions as specified in imagedesc.

	h : {nd_array}
	    the volume 1D amplitude histogram

	ih : {nd_array}
	    the volume 1D amplitude integrated histogram

	index : {int}
	    the amplitude (histogram index) where the integrated histogram
	    crosses the 'threshold' provided.

	Examples
	--------

	>>> import numpy as NP
	>>> import _registration as reg
	>>> anat_desc = reg.load_anatMRI_desc()
	>>> image_anat, h, ih, index = reg.load_volume(anat_desc, imagename='ANAT1_V0001.img', debug=1)
	>>> index
	210


	"""

	# load MRI or fMRI volume and return an autoscaled 8 bit image.
	# autoscale is using integrated histogram to deal with outlier high amplitude voxels
	if imagename == None:
	    # imagename of none means to create a blank image
    	    ImageVolume = NP.zeros(imagedesc['layers']*imagedesc['rows']*imagedesc['cols'],
			    dtype=NP.uint16).reshape(imagedesc['layers'], imagedesc['rows'], imagedesc['cols'])
	else:
    	    ImageVolume = NP.fromfile(imagename,
			    dtype=NP.uint16).reshape(imagedesc['layers'], imagedesc['rows'], imagedesc['cols']);

	# the mat (voxel to physical) matrix
	M = NP.eye(4, dtype=NP.float64);
	# for now just the sample size (mm units) in x, y and z
	M[0][0] = imagedesc['sample_x']
	M[1][1] = imagedesc['sample_y']
	M[2][2] = imagedesc['sample_z']
	# dimensions 
	D = NP.zeros(3, dtype=NP.int32);
	# Gaussian kernel - fill in with build_fwhm() 
	F = NP.zeros(3, dtype=NP.float64);
	D[0] = imagedesc['rows']
	D[1] = imagedesc['cols']
	D[2] = imagedesc['layers']

	if imagename == None:
	    # no voxels to scale to 8 bits
    	    ImageVolume = ImageVolume.astype(NP.uint8)
	    image = {'data' : ImageVolume, 'mat' : M, 'dim' : D, 'fwhm' : F}
    	    return image

	# 8 bit scale with threshold clip of the volume integrated histogram
	max = ImageVolume.max()
	min = ImageVolume.min()
	ih  = NP.zeros(max-min+1, dtype=NP.float64);
	h   = NP.zeros(max-min+1, dtype=NP.float64);
	if threshold <= 0:
	    threshold = 0.999
	elif threshold > 1.0:
	    threshold = 1.0
	# get the integrated histogram of the volume and get max from 
	# the threshold crossing in the integrated histogram 
	index  = R.register_image_threshold(ImageVolume, h, ih, threshold)
	scale  = 255.0 / (index-min)
	# generate the scaled 8 bit image
	images = (scale*(ImageVolume.astype(NP.float)-min))
	images[images>255] = 255 
	image = {'data' : images.astype(NP.uint8), 'mat' : M, 'dim' : D, 'fwhm' : F}
	if debug == 1:
    	    return image, h, ih, index
        else:
    	    return image



#
#  ---- demo/debug routines  ---- 
#

def load_anatMRI_desc():
	# this is for demo on the test MRI and fMRI volumes
	rows   = 256
	cols   = 256
	layers = 90
	xsamp  = 0.9375
	ysamp  = 0.9375
	zsamp  = 1.5
	desc = {'rows' : rows, 'cols' : cols, 'layers' : layers, 
		'sample_x' : xsamp, 'sample_y' : ysamp, 'sample_z' : zsamp}
	return desc

def load_fMRI_desc():
	# this is for demo on the test MRI and fMRI volumes
	rows   = 64
	cols   = 64
	layers = 28
	xsamp  = 3.75
	ysamp  = 3.75
	zsamp  = 5.0
	desc = {'rows' : rows, 'cols' : cols, 'layers' : layers, 
		'sample_x' : xsamp, 'sample_y' : ysamp, 'sample_z' : zsamp}
	return desc

def read_fMRI_directory(path):
	files_fMRI = glob.glob(path)
	return files_fMRI


def check_alignment(image1, image2, imdata, method='ncc', lite=0, smhist=0, 
		    alpha=0.0, beta=0.0, gamma=0.0, Tx=0, Ty=0, Tz=0, ret_histo=0):
	            
	#
	# to test the cost function and view the joint histogram
	# for 2 images. used for debug
	#
	imdata['parms'][0] = alpha
	imdata['parms'][1] = beta
	imdata['parms'][2] = gamma
	imdata['parms'][3] = Tx
	imdata['parms'][4] = Ty
	imdata['parms'][5] = Tz
	M = build_rotate_matrix(imdata['parms'])
	optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)

	if ret_histo:
	    cost, joint_histogram = optimize_function(imdata['parms'], optfunc_args)
	    return cost, joint_histogram 
    	else:
	    cost = optimize_function(imdata['parms'], optfunc_args)
	    return cost

def build_scale_image(image, scale):
	#
	# rescale the 'mat' (voxel to physical mapping matrix) 
	#
	(layers, rows, cols) = image['data'].shape
	M = image['mat'] * scale
	# dimensions 
	D = NP.zeros(3, dtype=NP.int32);
	# Gaussian kernel - fill in with build_fwhm() 
	F = NP.zeros(3, dtype=NP.float64);
	Z = NP.zeros(3, dtype=NP.float64);
	D[0] = rows/scale
	D[1] = cols/scale
	D[2] = layers/scale
    	image2 = NP.zeros(D[2]*D[1]*D[0], dtype=NP.uint8).reshape(D[2], D[0], D[1]);
	mode = 1;
	R.register_volume_resample(image['data'], image2, Z, scale, mode)
	scaled_image = {'data' : image2, 'mat' : M, 'dim' : D, 'fwhm' : F}
	return scaled_image


def demo_MRI_volume_align(scale=2, alpha=3.0, beta=4.0, gamma=5.0, Tx = 0.0, Ty = 0.0, Tz = 0.0):
	"""
	demo with (must have file ANAT1_V0001.img)

	image1, image2, imdata = reg.demo_MRI_volume_align()
	x = reg.python_coreg(image1, image2, imdata, method='ncc', lite=1) 
	image2r = reg.remap_image(image2, x, resample='cubic')
	image2rz = reg.resize_image(image2r, image1['mat'])


	slice1 = image1['data'][45, :, :]
	slice2 = image2['data'][45/2, :, :]
	slice2r = image2r['data'][45/2, :, :]
	slice2rz = image2rz['data'][45, :, :]

	pylab.figure(1)
	pylab.bone()
	pylab.imshow(slice1)
	pylab.imshow(slice1)
	pylab.figure(2)
	pylab.imshow(slice2)
	pylab.figure(3)
	pylab.imshow(slice2r)
	pylab.figure(4)
	pylab.imshow(slice2rz)
	pylab.show()

	"""
	#
	# this is for coreg MRI / fMRI scale test. The volume is anatomical MRI.
	# the image is rotated in 3D. after rotation the image is scaled.  
	#

	anat_desc = load_anatMRI_desc()
	image1 = load_volume(anat_desc, imagename='ANAT1_V0001.img')
	image2 = load_volume(anat_desc, imagename=None)
	imdata = build_structs()
	image1['fwhm'] = build_fwhm(image1['mat'], imdata['step'])
	image2['fwhm'] = build_fwhm(image2['mat'], imdata['step'])
	imdata['parms'][0] = alpha
	imdata['parms'][1] = beta
	imdata['parms'][2] = gamma
	imdata['parms'][3] = Tx
	imdata['parms'][4] = Ty
	imdata['parms'][5] = Tz
	M = build_rotate_matrix(imdata['parms'])
	# rotate volume. linear interpolation means the volume is low pass filtered
	R.register_linear_resample(image1['data'], image2['data'], M, imdata['step'])
	# subsample volume
	image3 = build_scale_image(image2, scale)
	return image1, image3, imdata

def demo_rotate_fMRI_volume(fMRIVol, x): 
	#
	# return rotated fMRIVol. the fMRIVol is already loaded, and gets rotated
	#

	desc = load_fMRI_desc()
	image = load_volume(desc, imagename=None)
	imdata = build_structs()
	image['fwhm'] = build_fwhm(image['mat'], imdata['step'])
	imdata['parms'][0] = x[0]  # alpha
	imdata['parms'][1] = x[1]  # beta
	imdata['parms'][2] = x[2]  # gamma
	imdata['parms'][3] = x[3]  # Tx
	imdata['parms'][4] = x[4]  # Ty
	imdata['parms'][5] = x[5]  # Tz
	M = build_rotate_matrix(imdata['parms'])
	# rotate volume. cubic spline interpolation means the volume is NOT low pass filtered
	R.register_cubic_resample(fMRIVol['data'], image['data'], M, imdata['step'])
	return image

def demo_MRI_coregistration(optimizer_method='powell', histo_method=1, smooth_histo=0, smooth_image=0, ftype=1):
	"""
	demo with (must have file ANAT1_V0001.img and fMRI directory fMRIData)

	measures, imageF_anat, fmri_series = reg.demo_MRI_coregistration()

	show results with

	In [59]: measures[25]['cost']
	Out[59]: -0.48607185

	In [60]: measures[25]['align_cost']
	Out[60]: -0.99514639

	In [61]: measures[25]['align_rotate']
	Out[61]:
	array([ 1.94480181,  5.64703989,  5.35002136, -5.00544405, -2.2712214, -1.42249691], dtype=float32)

	In [62]: measures[25]['rotate']
	Out[62]:
	array([ 1.36566341,  4.70644331,  4.68198586, -4.32256889, -2.47607017, -2.39173937], dtype=float32)


	"""

	# demo of alignment of fMRI series with anatomical MRI
	# in this demo, each fMRI volume is first perturbed (rotated, translated) 
	# by a random value. The initial registration is measured, then the optimal
	# alignment is computed and the registration measure made following the volume remap.
	# The fMRI registration is done with the first fMRI volume using normalized cross-correlation.
	# Each fMRI volume is rotated to the fMRI-0 volume and the series is ensemble averaged.
	# The ensemble averaged is then registered with the anatomical MRI volume using normalized mutual information.
	# The fMRI series is then rotated with this parameter. The alignments are done with 3D cubic splines.

	# read the anatomical MRI volume
	anat_desc = load_anatMRI_desc()
	imageF_anat = load_volume(anat_desc, imagename='ANAT1_V0001.img')
	# the sampling structure
	imdata = build_structs()
	# the volume filter
	imageF_anat['fwhm'] = build_fwhm(imageF_anat['mat'], imdata['step'])

	# read in the file list of the fMRI data
	metric_test = NP.dtype([('cost', 'f'),
                    	       ('align_cost', 'f'),
                    	       ('rotate', 'f', 6),
                    	       ('align_rotate', 'f', 6)])

	fMRIdata = read_fMRI_directory('fMRIData\*.img')
	fmri_desc = load_fMRI_desc()
	fmri_series = {}
	ave_fMRI_volume = NP.zeros(fmri_desc['layers']*fmri_desc['rows']*fmri_desc['cols'],
			  dtype=NP.float64).reshape(fmri_desc['layers'], fmri_desc['rows'], fmri_desc['cols'])
	count = 0
	number_volumes = len(fMRIdata)
	measures = NP.zeros(number_volumes, dtype=metric_test)
	# load and perturb (rotation, translation) the fMRI volumes
	for i in fMRIdata:
	    image = load_volume(fmri_desc, i)
	    # random perturbation of angle, translation for each volume beyond the first
	    if count == 0:
		image['fwhm'] = build_fwhm(image['mat'], imdata['step'])
	        fmri_series[count] = image
	        count = count + 1
	    else:
	        x = NP.random.random(6) - 0.5
	        x = 10.0 * x
	        fmri_series[count] = demo_rotate_fMRI_volume(image, x)
		measures[count]['rotate'][0:6] = x[0:6]
	        count = count + 1


	# load and register the fMRI volumes with volume_0 using normalized cross correlation metric
	imageF = fmri_series[0]
	if smooth_image:
	    image_F_xyz = filter_image_3D(imageF['data'], imageF['fwhm'], ftype)
	    imageF['data'] = image_F_xyz
	for i in range(1, number_volumes):
	    imageG = fmri_series[i]
	    # the measure prior to alignment 
	    measures[i]['cost'] = check_alignment(imageF, imageG, imdata, method='ncc',
			                          lite=histo_method, smhist=smooth_histo)
	    x = python_coreg(imageF, imageG, imdata, lite=histo_method, method='ncc',
			     opt_method=optimizer_method, smhist=smooth_histo, smimage=smooth_image)
	    measures[i]['align_rotate'][0:6] = x[0:6]
	    measures[i]['align_cost'] = check_alignment(imageF, imageG, imdata, method='ncc', 
		                             lite=histo_method, smhist=smooth_histo,
					     alpha=x[0], beta=x[1], gamma=x[2], Tx=x[3], Ty=x[4], Tz=x[5])


	# align the volumes and average them for co-registration with the anatomical MRI 
	ave_fMRI_volume = fmri_series[0]['data'].astype(NP.float64)
	for i in range(1, number_volumes):
	    image = fmri_series[i]
	    x[0:6] = measures[i]['align_rotate'][0:6]
	    # overwrite the fMRI volume with the aligned volume
	    fmri_series[i] = remap_image(image, x, resample='cubic')
	    ave_fMRI_volume = ave_fMRI_volume + fmri_series[i]['data'].astype(NP.float64)

	ave_fMRI_volume = (ave_fMRI_volume / float(number_volumes)).astype(NP.uint8)
	ave_fMRI_volume = {'data' : ave_fMRI_volume, 'mat' : imageF['mat'], 
			   'dim' : imageF['dim'], 'fwhm' : imageF['fwhm']}
	# register (using normalized mutual information) with the anatomical MRI
	if smooth_image:
	    image_F_anat_xyz = filter_image_3D(imageF_anat['data'], imageF_anat['fwhm'], ftype)
	    imageF_anat['data'] = image_F_anat_xyz
	x = python_coreg(imageF_anat, ave_fMRI_volume, imdata, lite=histo_method,
			 method='nmi', opt_method=optimizer_method, smhist=smooth_histo, smimage=smooth_image)
	print 'functional-anatomical align parameters '
	print x
	for i in range(number_volumes):
	    image = fmri_series[i]
	    # overwrite the fMRI volume with the anatomical-aligned volume
	    fmri_series[i] = remap_image(image, x, resample='cubic')

	return measures, imageF_anat, fmri_series


def demo_fMRI_resample(imageF_anat, fmri_series):
	resampled_fmri_series = {}
	number_volumes = len(fmri_series)
	for i in range(number_volumes):
	    resampled_fmri_series[i] = resize_image(fmri_series[i], imageF_anat['mat'])

	return resampled_fmri_series


