import math

import numpy as np
from scipy.special import erf
from scipy.ndimage import correlate1d
from scipy.optimize import fmin_powell, fmin_cg

import scipy.ndimage._register as reg

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

def resize_image(imageS, imageS_mat, imageR_mat):
    """
    zoom_image = resize_image(imageS, imageS_mat, imageR_mat)

    Fractional resample source_image to reference_image size. The
    resample is implemented with 3D cubic spline. The source
    imageS_mat is the 4x4 voxel-to-physical conversion matrix.
    
    Parameters 
    ----------

    imageS: {ndarray} 
        imageS is the source image to be resized.

    imageS_mat : {ndarray} 
        the 4x4 transform of the source image that maps voxel to physical.

    imageR_mat : {ndarray}
        the 4x4 transform of the destination image that maps voxel to physical.

    Returns 
    -------
    zoom_image : {ndarray}

    Examples
    --------

    >>> import _registration as reg
    >>> measures, image_anat, image_anat_mat, image_fmri_mat, fmri_series = reg.demo_MRI_coregistration()

    >>> resampled_fmri = reg.resize_image(fmri_series[10], image_fmri_mat, image_anat_mat)

    image 10 in the fmri_series is resampled from image_fmri_mat to image_anat coordinates

    """

    # get the zoom
    Z = imageS_mat.diagonal() / imageR_mat.diagonal()

    # new volume dimensions (rounded). D, imageS and Z are 3D and this is a vector element product
    D = (imageS.shape * Z + 0.5).astype(np.int16)

    # for the test data, set the xyz voxel sizes for fMRI volume. M is a 4x4 matrix.
    M = np.diag(imageS_mat.diagonal() / Z)    

    image = np.zeros((D[2],D[1],D[0]),np.uint8)
    
    mode  = 2
    scale = 0
    reg.register_volume_resample(imageS, image, Z, scale, mode)

    return image, M

def remap_image(image, parm_vector, resample='linear'):
    """
    remaped_image = remap_image(image, parm_vector, resample='linear')

    rotates and translates image using the 3 angles and 3 translations in the 6-dim
    parm_vector. The mapping is stored in the 4x4 M_inverse matrix from the get_inverse_mapping
    method.

    Parameters 
    ----------
    image : {ndarray} 
        image is the source image to be remapped. 

    parm_vector : {ndarray}
        parm_vector is the 6-dimensional vector (3 angles, 3 translations)
        generated from the rigid body registration.

    resample : {'linear', 'cubic'}, optional


    Returns 
    -------
    remaped_image : {ndarray}

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

    # allocate the zero image
    remaped_image = np.zeros(image.shape, dtype=np.uint8)

    step = np.array([1, 1, 1], dtype=np.int32)

    if resample == 'linear':
        # trilinear interpolation mapping.
        reg.register_linear_resample(image, remaped_image, M_inverse, step)
    elif resample == 'cubic':
        # tricubic convolve interpolation mapping. 
        reg.register_cubic_resample(image, remaped_image, M_inverse, step)

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
    >>> array = np.zeros(6, dtype=float)
    >>> M = reg.get_inverse_mappings(array)
    >>> M 

    array([
    [ 1.,  0.,  0.,  0.],
    [ 0.,  1.,  0.,  0.],
    [ 0.,  0.,  1.,  0.],
    [ 0.,  0.,  0.,  1.]])

    """
    # get the inverse mapping to rotate the G matrix to F space following registration
    # -parm_vector is the inverse angles and translations
    M_inverse = build_rotate_matrix(-parm_vector)
    return M_inverse

def register(image1, image1_mat, image2, image2_mat, multires=[4, 2], histo_fwhm=3, 
             ftype=1, lite=0, smhist=0, method='nmi', opt_method='hybrid',
	     optimize_function=None):

    """
    parm_vector = register(image1, image1_mat, image2, image2_mat, multires=[4, 2], histo_fwhm=3,
                             ftype=1, lite=0, smhist=0, method='nmi', opt_method='powell'):

    alignment of the two images (measured by mutual information or cross correlation) 
    using optimization search of 3 angle and 3 translation parameters. The optimization 
    uses either the Powell or Conjugate Gradient methods in the scipy optimization 
    package. The optimal rigid body parameter is returned.

    Parameters 
    ----------
    image1 : {nd_array} 
        image1 is the source image to be remapped during the registration. 
    image1_mat : {nd_array} 
        image1_mat is the source image MAT 
    image2 : {nd_array} 
        image2 is the reference image that image1 gets mapped to. 
    image2_mat : {nd_array} 
        image2_mat is the source image MAT 
    multires: {list}, optional
        the volume subsample values for each pass of the registration.
	the default is 2 passes with subsample 4 in pass 1 and subsample 2 in pass 2
    histo_fwhm : {int}, optional
        used for the filter kernel in the low pass filter of the joint histogram 
    ftype : {0, 1}, optional
        flag for type of low pass filter. 0 is Gauss-Spline
        1 is pure Gauss. Sigma determined from volume sampling info.
    lite : {0, 1}, optional
        lite of 1 is to jitter both images during resampling. 0
        is to not jitter. jittering is for non-aliased volumes.
    smhist: {0, 1}, optional
        flag for joint histogram low pass filtering. 0 for no filter,
        1 for do filter.
    method: {'nmi', 'mi', 'ncc', 'ecc', 'mse'}, optional
        flag for type of registration metric. nmi is normalized mutual
        information; mi is mutual information; ecc is entropy cross
        correlation; ncc is normalized cross correlation. mse is mean
	squared error.
    opt_method: {'powell', 'cg', 'hybrid'}, optional
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

    >>> image1, image2, fwhm, improc = reg.demo_build_dual_volumes()
    >>> parm_vector = register(image1, image2, fwhm, improc)

    """

    # do the parameter validity checking. this is specific to this 3D registration.
    # make sure the image is 3D and the mats are 4x4 with nonzero diagonal

    if image1.ndim != 3:
        raise ValueError, "Image 1 is not 3 dimensional"

    if image2.ndim != 3:
        raise ValueError, "Image 2 is not 3 dimensional"

    if image1.dtype != np.uint8:
        raise ValueError, "Image 1 is not 8 bit (required for joint histogram)"

    if image2.dtype != np.uint8:
        raise ValueError, "Image 2 is not 8 bit (required for joint histogram)"

    if image1_mat.shape != (4,4):
        raise ValueError, "Image1 MAT is not 4x4"

    if image2_mat.shape != (4,4):
        raise ValueError, "Image2 MAT is not 4x4"

    if (np.diag(image1_mat)).prod() == 0:
        raise ValueError, "Image1 MAT has a 0 on the diagonal"

    if (np.diag(image2_mat)).prod() == 0:
        raise ValueError, "Image2 MAT has a 0 on the diagonal"

    if opt_method=='hybrid' and np.size(multires) != 2:
        raise ValueError, "hybrid method must be 2 pass registration"

    if ftype != 0 and ftype != 1: 
        raise ValueError, "choose filter type 0 or 1 only"

    if lite != 0 and lite != 1: 
        raise ValueError, "choose histogram generation type 0 or 1 only"

    if smhist != 0 and smhist != 1: 
        raise ValueError, "choose histogram smoothing type 0 or 1 only"

    if method != 'nmi' and method != 'mi'  and method != 'ncc'\
                       and method != 'ecc' and method != 'mse':
        raise ValueError, "choose cost method nmi, mi, ecc, mse, ncc"

    if opt_method != 'powell' and opt_method != 'cg'  and opt_method != 'hybrid':
        raise ValueError, "only optimize methods powell, cg or hybrid are supported"

    # default is to use the cost_function I provided.
    # this shows you can override this but the parameters will have to
    # be changed for the new cost function if it is different

    if optimize_function == None:
        optimize_function = cost_function

    parm_vector = multires_registration(optimize_function, image1, image1_mat, image2, image2_mat,
		                        multires, histo_fwhm, lite, smhist, method, opt_method)

    return parm_vector

def multires_registration(optimize_function, image1, image1_mat, image2, image2_mat,
		          multires, histo_fwhm, lite, smhist, method, opt_method):

    """

    to be called by register() which does parameter validation 

    Parameters 
    ----------
    image1 : {nd_array} 
        image1 is the source image to be remapped during the registration. 
    image1_mat : {nd_array} 
        image1_mat is the source image MAT 
    image2 : {nd_array} 
        image2 is the reference image that image1 gets mapped to. 
    image2_mat : {nd_array} 
        image2_mat is the source image MAT 
    multires: {list}, optional
        the volume subsample values for each pass of the registration.
	the default is 2 passes with subsample 4 in pass 1 and subsample 2 in pass 2
    histo_fwhm : {int}, optional
        used for the filter kernel in the low pass filter of the joint histogram 
    ftype : {0, 1}, optional
        flag for type of low pass filter. 0 is Gauss-Spline
        1 is pure Gauss. Sigma determined from volume sampling info.
    lite : {0, 1}, optional
        lite of 1 is to jitter both images during resampling. 0
        is to not jitter. jittering is for non-aliased volumes.
    smhist: {0, 1}, optional
        flag for joint histogram low pass filtering. 0 for no filter,
        1 for do filter.
    method: {'nmi', 'mi', 'ncc', 'ecc', 'mse'}, optional
        flag for type of registration metric. nmi is normalized mutual
        information; mi is mutual information; ecc is entropy cross
        correlation; ncc is normalized cross correlation. mse is mean
	squared error.
    opt_method: {'powell', 'cg', 'hybrid'}, optional
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

    (calling this from register which optionally filters image2)
    >>> import numpy as NP
    >>> import _registration as reg
    >>> image1, mat1, image2, mat2 = reg.demo_build_dual_volumes()
    >>> parm_vector = register(image1, image2, imdata)

    """
    ret_histo=0
    step = np.array([1, 1, 1], dtype=np.int32)
    fwhm = np.zeros(2, dtype=np.int32)
    # make the step a scalar to can put in a multi-res loop
    loop = range(np.size(multires))
    # 6-D zero vector
    x = np.zeros(6, dtype=np.float64);
    # the kernel fwhm value for the x and y joint histogram filter
    fwhm[:] = histo_fwhm
    for i in loop:
	# this is the volume subsample
	step[:] = multires[i]
	# optfunc_args is specific to the cost_function in this file
	# this will need to change if you use another optimize_function.
        optfunc_args = (image1, image1_mat, image2, image2_mat, step, histo_fwhm,
			lite, smhist, method, ret_histo)
        p_args = (optfunc_args,)
        if opt_method=='powell':
            print 'POWELL multi-res registration step size ', step
            print 'vector ', x
            x = fmin_powell(optimize_function, x, args=p_args, callback=callback_powell)
        elif opt_method=='cg':
            print 'CG multi-res registration step size ', step
            print 'vector ', x
            x = fmin_cg(optimize_function, x, args=p_args, callback=callback_cg) 
        elif opt_method=='hybrid':
            if i==0:
                print 'Hybrid POWELL multi-res registration step size ', step
                print 'vector ', x
                lite = 0
                optfunc_args = (image1, image1_mat, image2, image2_mat, step, histo_fwhm,
                                lite, smhist, method, ret_histo)
                p_args = (optfunc_args,)
                x = fmin_powell(optimize_function, x, args=p_args, callback=callback_powell) 
            elif i==1:
                print 'Hybrid CG multi-res registration step size ', step
                print 'vector ', x
                lite = 1
                optfunc_args = (image1, image1_mat, image2, image2_mat, step, histo_fwhm,
                                lite, smhist, method, ret_histo)
                p_args = (optfunc_args,)
                x = fmin_cg(optimize_function, x, args=p_args, callback=callback_cg) 

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

def smooth_kernel(fwhm, x, pixel_scale=8.0, ktype=1):
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
    >>> p = np.ceil(2*fwhm).astype(int)
    >>> x = np.array(range(-p, p+1))
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
    s   = np.square((fwhm/math.sqrt(pixel_scale*math.log(2.0)))) + eps
    if ktype==1:
        # from SPM: Gauss kernel convolved with 1st degree B spline
        w1 = 0.5 * math.sqrt(2.0/s)
        w2 = -0.5 / s
        w3 = math.sqrt((s*math.pi) /2.0)
        kernel = 0.5*(erf(w1*(x+1))*(x+1) + erf(w1*(x-1))*(x-1)
                      - 2.0*erf(w1*(x))*(x) + w3*(np.exp(w2*np.square(x+1))) 
                      + np.exp(w2*(np.square(x-1)))
                      - 2.0*np.exp(w2*np.square(x)))
        kernel[kernel<0] = 0
        kernel = kernel / kernel.sum()  
    else:
        # Gauss kernel 
        kernel = (1.0/math.sqrt(2.0*math.pi*s)) * np.exp(-np.square(x)/(2.0*s)) 
        kernel = kernel / kernel.sum()  

    return kernel

def filter_image_3D(imageRaw, fwhm, ftype=2, give_2D=0):
    """
    image_F_xyz = filter_image_3D(imageRaw, fwhm, ftype=2):
    does 3D separable digital filtering using scipy.ndimage.correlate1d

    Parameters 
    ----------
    imageRaw : {nd_array}
        the unfiltered 3D volume image
    fwhm : {int}
        used for kernel width. this is 3 elements (one for each dimension)
    ktype: {1, 2}, optional
        kernel type. 1 is Gauss convoled with spline (SPM), 2 is Gauss

    Returns 
    -------
    image_F_xyz : {nd_array}
        3D filtered volume image

    Examples
    --------

    >>> import _registration as reg
    >>> image1, image2, imdata = reg.demo_build_dual_volumes()
    >>> ftype = 1
    >>> image_Filter_xyz = filter_image_3D(image, fwhm, ftype)
    >>> image1['data'] = image_Filter_xyz
    """

    p = np.ceil(2*fwhm).astype(int)
    x = np.array(range(-p[0], p[0]+1))
    kernel_x = smooth_kernel(fwhm[0], x, ktype=ftype)

    x = np.array(range(-p[1], p[1]+1))
    kernel_y = smooth_kernel(fwhm[1], x, ktype=ftype)

    x = np.array(range(-p[2], p[2]+1))
    kernel_z = smooth_kernel(fwhm[2], x, ktype=ftype)

    output=None
    # 3D filter in 3 1D separable stages. keep the image
    # names at each stage separate in case you need them
    # for example may need an image that is 2D slice filtered only
    axis = 0
    image_F_x   = correlate1d(imageRaw,   kernel_x, axis, output)
    axis = 1
    image_F_xy  = correlate1d(image_F_x,  kernel_y, axis, output)
    axis = 2
    image_F_xyz = correlate1d(image_F_xy, kernel_z, axis, output)

    if give_2D==0:
        return image_F_xyz  
    else:
        return image_F_xyz, image_F_xy


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
    >>> image1['fwhm'] = reg.build_fwhm(image1['mat'], imdata['step'])

    """
    # M contains the voxel to physical mapping
    view_3x3 = np.square(M[0:3, 0:3])
    # sum the elements in the first row
    vxg = np.sqrt(view_3x3.sum(axis=0))
    # assumes that voxel sampling is the same for xyz as S is the step
    size = np.array([1,1,1])*S[0]
    x = np.square(size) - np.square(vxg)
    # clip
    x[x<0] = 0
    fwhm = np.sqrt(x) / vxg
    # pathology when stepsize = 1 for MAT equal to the identity matrix
    fwhm[fwhm==0] = 1
    # return the 3D Gaussian kernel width (xyz)
    return fwhm 

def cost_function(x, optfunc_args):
    """
    cost = cost_function(x, optfunc_args)    --- OR ---
    cost, joint_histogram = cost_function(x, optfunc_args)   

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
    >>> image1['fwhm'] = reg.build_fwhm(image1['mat'], imdata['step'])
    >>> image2['fwhm'] = reg.build_fwhm(image2['mat'], imdata['step'])
    >>> method = 'ncc'
    >>> lite = 1
    >>> smhist = 0
    >>> ret_histo = 1
    >>> optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)
    >>> x = np.zeros(6, dtype=np.float64)
    >>> return cost, joint_histogram = reg.cost_function(x, optfunc_args)


    """

    image_F       = optfunc_args[0]
    image_F_mat   = optfunc_args[1]
    image_G       = optfunc_args[2]
    image_G_mat   = optfunc_args[3]
    sample_vector = optfunc_args[4]
    fwhm          = optfunc_args[5]
    do_lite       = optfunc_args[6]
    smooth        = optfunc_args[7]
    method        = optfunc_args[8]
    ret_histo     = optfunc_args[9]

    rot_matrix = build_rotate_matrix(x)
    cost = 0.0
    epsilon = 2.2e-16 
    # image_G is base image
    # image_F is the to-be-rotated image
    # rot_matrix is the 4x4 constructed (rigid body) transform matrix
    # sample_vector is the subsample vector for x-y-z

    F_inv = np.linalg.inv(image_F_mat)
    composite = np.dot(F_inv, image_G_mat)
    composite = np.dot(composite, rot_matrix)

    if method == 'mse':
        #
        # mean squard error method
        #

        # allocate the zero image
        #(layers, rows, cols) = image_F.shape
        remap_image_F = np.zeros(image_F.shape, dtype=np.uint8)
        # trilinear interpolation mapping.
        reg.register_linear_resample(image_F, remap_image_F, composite, sample_vector)
        cost = (np.square(image_G-remap_image_F)).mean()
	# cost is min when G and F are aligned so keep cost positive

        return cost

    else:
        #
        # histogram-based methods (nmi, ncc, mi, ecc)
        #

        # allocate memory for 2D histogram
        joint_histogram = np.zeros([256, 256], dtype=np.float64)

        if do_lite: 
            reg.register_histogram_lite(image_F, image_G, composite, sample_vector, joint_histogram)
        else:
            reg.register_histogram(image_F, image_G, composite, sample_vector, joint_histogram)

        # smooth the histogram
        if smooth: 
            p = np.ceil(2*fwhm).astype(int)
            x = np.array(range(-p, p+1))
            hkernel = smooth_kernel(fwhm, x)
            output=None
            # 2D filter in 1D separable stages using the same kernel. SPM
	    # has options for a 2D fwhm kernel yet only uses 1 element
            axis = 0
            joint_histogram = correlate1d(joint_histogram, hkernel, axis, output)
            axis = 1
            joint_histogram = correlate1d(joint_histogram, hkernel, axis, output)

        joint_histogram += epsilon # prevent log(0) 
        # normalize the joint histogram
        joint_histogram /= joint_histogram.sum() 
        # get the marginals
        marginal_col = joint_histogram.sum(axis=0)
        marginal_row = joint_histogram.sum(axis=1)

        if method == 'mi':
            # mutual information
            marginal_outer = np.outer(marginal_col, marginal_row)
            H = joint_histogram * np.log(joint_histogram / marginal_outer)  
            mutual_information = H.sum()
            cost = -mutual_information

        elif method == 'ecc':
            # entropy correlation coefficient 
            marginal_outer = np.outer(marginal_col, marginal_row)
            H = joint_histogram * np.log(joint_histogram / marginal_outer)  
            mutual_information = H.sum()
            row_entropy = marginal_row * np.log(marginal_row)
            col_entropy = marginal_col * np.log(marginal_col)
            ecc  = -2.0*mutual_information/(row_entropy.sum() + col_entropy.sum())
            cost = -ecc

        elif method == 'nmi':
            # normalized mutual information
            row_entropy = marginal_row * np.log(marginal_row)
            col_entropy = marginal_col * np.log(marginal_col)
            H = joint_histogram * np.log(joint_histogram)  
            nmi  = (row_entropy.sum() + col_entropy.sum()) / (H.sum())
            cost = -nmi

        elif method == 'ncc':
            # cross correlation from the joint histogram 
            r, c = joint_histogram.shape
            i = np.array(range(1,c+1))
            j = np.array(range(1,r+1))
            m1 = (marginal_row * i).sum()
            m2 = (marginal_col * j).sum()
            sig1 = np.sqrt((marginal_row*(np.square(i-m1))).sum())
            sig2 = np.sqrt((marginal_col*(np.square(j-m2))).sum())
            [a, b] = np.mgrid[1:c+1, 1:r+1]
            a = a - m1
            b = b - m2
            # element multiplies in the joint histogram and grids
            H = ((joint_histogram * a) * b).sum()
            ncc  = H / (np.dot(sig1, sig2)) 
            cost = -ncc

        if ret_histo:
            return cost, joint_histogram 
        else:
            return cost


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
    >>> x = np.zeros(6, dtype=np.float64)
    >>> M = reg.build_rotate_matrix(x)
    >>> M 
    array([[ 1.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.],
           [ 0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  1.]])


    """

    R1 = np.zeros([4,4], dtype=np.float64);
    R2 = np.zeros([4,4], dtype=np.float64);
    R3 = np.zeros([4,4], dtype=np.float64);
    T  = np.eye(4, dtype=np.float64);

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

    rot_matrix = np.dot(T, R1);
    rot_matrix = np.dot(rot_matrix, R2);
    rot_matrix = np.dot(rot_matrix, R3);

    return rot_matrix


def build_test_volume(imagedesc, S=[1500.0, 2500.0, 1000.0]):

    """
    build a 3D Gaussian volume. user passes in image dims in imagedesc
    the sigma for each axis is S[3] where 0=z, 1=y, 2=x

    volume3D = build_test_volume(imagedesc, S)

    Parameters 
    ----------
    imagedesc : {dictionary}
        volume dimensions and sampling

    S : {tuple}
        the Gaussian sigma for Z, Y and X

    Returns 
    -------

    volume3D : {nd_array}
        the 3D volume for testing

    """
    layers = imagedesc['layers']
    rows   = imagedesc['rows']
    cols   = imagedesc['cols']

    L = layers/2
    R = rows/2
    C = cols/2

    # build coordinates for 3D Gaussian volume
    # coordinates are centered at (0, 0, 0)
    [a, b, c] = np.mgrid[-L:L, -R:R, -C:C]

    sigma    = np.array([S[0], S[1], S[2]])
    aa       = (np.square(a))/sigma[0]
    bb       = (np.square(b))/sigma[1]
    cc       = (np.square(c))/sigma[2]
    volume3D = (255.0*np.exp(-(aa + bb + cc))).astype(np.uint8)

    return volume3D



def load_volume(imagedesc, imagename=None):

    """

    returns an nd_array from the filename or blank image. used for testing. 

    Parameters 
    ----------
    imagedesc : {dictionary} 
        imagedesc is the descriptor of the image to be read. 

    imagename : {string} : optional
        name of image file. No name creates a blank image that is used for creating
        a rotated test image or image rescaling.

    Returns 
    -------
    image : {nd_array}
        the volume data assoicated with the filename or a blank volume of the same
        dimensions as specified in imagedesc.

    M : {nd_array}
        the voxel-to-physical affine matrix (mat)

    Examples
    --------

    >>> import _registration as reg
    >>> anat_desc = reg.load_anatMRI_desc()
    >>> image, M = reg.load_volume(anat_desc, imagename='ANAT1_V0001.img')


    """

    # load MRI or fMRI volume and return an autoscaled 8 bit image.
    # autoscale is using integrated histogram to deal with outlier high amplitude voxels
    if imagename == None:
        # imagename of none means to create a blank image
        image = np.zeros([imagedesc['layers'],imagedesc['rows'],imagedesc['cols']],dtype=np.uint16)
    else:
        image = np.fromfile(imagename,
                   dtype=np.uint16).reshape(imagedesc['layers'], imagedesc['rows'], imagedesc['cols']);

    # the mat (voxel to physical) matrix
    M = np.eye(4, dtype=np.float64);
    # for now just the sample size (mm units) in x, y and z
    M[0][0] = imagedesc['sample_x']
    M[1][1] = imagedesc['sample_y']
    M[2][2] = imagedesc['sample_z']

    if imagename == None:
        # no voxels to scale to 8 bits
        image = image.astype(np.uint8)

    return image, M



def scale_image(image, max_amp=255, image_type=np.uint8, threshold=0.999, fetch_ih=0):

    """
    scale and threshold clip the volume using the integrated histogram
    to set the high threshold

    Parameters 
    ----------
    image : {nd_array}
        raw unscaled volume

    max_amp : int (default 255)
        the maximum value of the scaled image

    image_type : nd_array dtype (default uint8)
        the type of the volume to return.

    threshold : float (default 0.999)
        the value of the normalized integrated histogram
	that when reached sets the high threshold index

    Returns 
    -------
    image : {nd_array}
        the scaled volume
    ih : {nd_array}
        the integrated histogram. can be used for image display 
	purpose (histogram equalization)

    """

    max = image.max()
    min = image.min()
    if max == 0 and min == 0:
        raise ValueError, "Zero image. cannot be scaled"

    # need range of pixels for the number of bins
    h, edges = np.histogram(image, bins=(max-min))
    ih = (np.cumsum(h)).astype(np.float64)
    # normalize the integrated histogram
    ih = ih / ih.max()
    indices = np.where(ih >= threshold)
    # wind up getting all the indices where the ih >= threshold
    # and only need the first index. tuple has one nd_array and
    # get the 0 element from it ([0][0])
    index   = indices[0][0]
    scale   = float(max_amp) / (index-min)
    image   = (scale*(image.astype(np.float)-min))
    image[image>max_amp] = max_amp
    # down type. usually will go from float to 8 bit (needed for the 8 bit joint histogram)
    image = image.astype(image_type)

    if fetch_ih == 1:
        return image, ih
    else:
        return image


def check_alignment(image1, image1_mat, image2, image2_mat, histo_fwhm=3, method='ncc', lite=0,
                    smhist=0, alpha=0.0, beta=0.0, gamma=0.0, Tx=0, Ty=0, Tz=0, ret_histo=0):
                    
    """
    test the cost function and (optional) view the joint histogram. can be used
    during intra-modal registration to measure the current alignment (return
    the cross correlation). would measure before and after registration



    """

    # do the parameter validity checking. this is specific to this 3D registration.
    # make sure the image is 3D and the mats are 4x4 with nonzero diagonal

    if image1.ndim != 3:
        raise ValueError, "Image 1 is not 3 dimensional"

    if image2.ndim != 3:
        raise ValueError, "Image 2 is not 3 dimensional"

    if image1.dtype != np.uint8:
        raise ValueError, "Image 1 is not 8 bit (required for joint histogram)"

    if image2.dtype != np.uint8:
        raise ValueError, "Image 2 is not 8 bit (required for joint histogram)"

    if image1_mat.shape != (4,4):
        raise ValueError, "Image1 MAT is not 4x4"

    if image2_mat.shape != (4,4):
        raise ValueError, "Image2 MAT is not 4x4"

    if (np.diag(image1_mat)).prod() == 0:
        raise ValueError, "Image1 MAT has a 0 on the diagonal"

    if (np.diag(image2_mat)).prod() == 0:
        raise ValueError, "Image2 MAT has a 0 on the diagonal"

    if method != 'nmi' and method != 'mi'  and method != 'ncc'\
                       and method != 'ecc' and method != 'mse':
        raise ValueError, "choose cost method nmi, mi, ecc, mse, ncc"

    P    = np.zeros(6, dtype=np.float64);
    P[0] = alpha
    P[1] = beta
    P[2] = gamma
    P[3] = Tx
    P[4] = Ty
    P[5] = Tz

    step = np.array([1, 1, 1], dtype=np.int32)
    optfunc_args = (image1, image1_mat, image2, image2_mat, step, histo_fwhm, lite,
		    smhist, method, ret_histo)
			
    if ret_histo:
        cost, joint_histogram = cost_function(P, optfunc_args)
        return cost, joint_histogram 
    else:
        cost = cost_function(P, optfunc_args)
        return cost



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



def build_scale_volume(image, mat, scale):
    #
    # rescale the 'mat' (voxel to physical mapping matrix) 
    #
    M = mat * scale
    (layers, rows, cols) = image.shape
    # dimensions 
    D = np.zeros(3, dtype=np.int32);
    Z = np.zeros(3, dtype=np.float64);
    D[0] = rows/scale
    D[1] = cols/scale
    D[2] = layers/scale
    image2 = np.zeros([D[2], D[0], D[1]], dtype=np.uint8)
    mode = 1;
    reg.register_volume_resample(image, image2, Z, scale, mode)
    return image2, M


def demo_build_dual_volumes(scale=2, alpha=3.0, beta=4.0, gamma=5.0, Tx = 0.0, Ty = 0.0, Tz = 0.0):
    """
    demo with (must have file ANAT1_V0001.img)
    builds a volume and a scaled-rotated version for coreg testing 

    image1, mat1, image2, mat2 = reg.demo_build_dual_volumes()
    x = reg.register(image1, mat1, image2, mat2, method='ncc', lite=1) 
    image2r = reg.remap_image(image2, x, resample='cubic')
    image2rz = reg.resize_image(image2r, mat1)

    """
    #
    # this is for coreg MRI / fMRI scale test. The volume is anatomical MRI.
    # the image is rotated in 3D. after rotation the image is scaled.  
    #

    step = np.array([1, 1, 1], dtype=np.int32)
    anat_desc = load_anatMRI_desc()
    image1, mat1 = load_volume(anat_desc, imagename='ANAT1_V0001.img')
    image2, mat2 = load_volume(anat_desc, imagename=None)
    image1 = scale_image(image1)
    P    = np.zeros(6, dtype=np.float64);
    P[0] = alpha
    P[1] = beta
    P[2] = gamma
    P[3] = Tx
    P[4] = Ty
    P[5] = Tz
    M = build_rotate_matrix(P)
    # rotate volume. linear interpolation means the volume is low pass filtered
    reg.register_linear_resample(image1, image2, M, step)
    # subsample volume
    image2, mat2 = build_scale_volume(image2, mat2, scale)
    return image1, mat1, image2, mat2

def tests(image1, mat1, image2, mat2): 

    # for same image, using the lite method the off-diagonal is zero
    cost, joint = reg.check_alignment(image1, mat1, image2, mat2, ret_histo=1, lite=1)
    my_diag = joint.diagonal()
    Z = np.diag(my_diag)
    W = joint - Z
    W[abs(W) < 1e-10] = 0.0

    if W.max() != 0.0 and W.min() != 0.0:
	print 'joint histogram is not diagonal '
    if abs(cost) < 0.99: 
	print 'cross correlation is too small'

    # for same image, not using the lite method the off-diagonal is non-zero 
    cost, joint = reg.check_alignment(image1, mat1, image2, mat2, ret_histo=1, lite=0)
    my_diag = joint.diagonal()
    Z = np.diag(my_diag)
    W = joint - Z
    W[abs(W) < 1e-10] = 0.0

    if W.max() == 0.0 and W.min() == 0.0:
	print 'joint histogram is diagonal and needs off-diagonals'
    if abs(cost) < 0.99: 
	print 'cross correlation is too small'

    # call w/o returning the joint histogram 
    cost = reg.check_alignment(image1, mat1, image2, mat2, ret_histo=0, lite=1)
    if abs(cost) < 0.99: 
	print 'cross correlation is too small'

    cost = reg.check_alignment(image1, mat1, image2, mat2, ret_histo=0, lite=0)
    if abs(cost) < 0.99: 
	print 'cross correlation is too small'


    image1 = np.zeros([64,64,64],np.uint8)
    image2 = np.zeros([64,64,64],np.uint8)
    image3 = np.zeros([64,64],np.uint8)
    mat1   = np.eye(4)
    mat2   = np.eye(3)
    mat3   = np.zeros([4,4])
    # test with wrong dim image, wrong dim mat and mat with zeros on diagonal

    # wrong image dim
    assertRaises(ValueError, check_alignment, image1, mat1, image3, mat1)
    # wrong mat dim
    assertRaises(ValueError, check_alignment, image1, mat1, image2, mat2)
    # mat with zeros on diagonal
    assertRaises(ValueError, check_alignment, image1, mat1, image2, mat3)












def demo_rotate_fMRI_volume(fMRI_volume, desc, x): 
    #
    # return rotated fMRIVol.
    #

    image = load_volume(desc, imagename=None)
    image = scale_image(image)
    step = np.array([1, 1, 1], dtype=np.int32)
    M = build_rotate_matrix(x)
    # rotate volume. cubic spline interpolation means the volume is NOT low pass filtered
    reg.register_cubic_resample(fMRI_volume, image, M, step)

    return image

def demo_MRI_coregistration(anatfile, funclist, optimizer_method='powell', 
                            histo_method=1, smooth_histo=0, smooth_image=0, 
                            ftype=1):
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
    imageF_anat, anat_mat = load_volume(anat_desc, imagename=anatfile)
    imageF = imageF_anat.copy()
    # the sampling structure
    step = np.array([1, 1, 1], dtype=np.int32)
    # the volume filter
    imageF_anat_fwhm = build_fwhm(anat_mat, step)


    # allocate the structure for the processed fMRI array
    metric_test = np.dtype([('cost', 'f'),
                            ('align_cost', 'f'),
                            ('rotate', 'f', 6),
                            ('align_rotate', 'f', 6)])
    # allocate the empty dictionary that will contain metrics and aligned volumes 
    fmri_series = {}

    # read in the file list of the fMRI data
    fMRIdata = read_fMRI_directory('fMRIData\*.img')
    fmri_desc = load_fMRI_desc()
    image_fmri, fmri_mat = load_volume(fmri_desc, fMRIdata[0])

    # one time build of the fwhm that is used to build the filter kernels
    anat_fwhm = build_fwhm(anat_mat, step)
    fmri_fwhm = build_fwhm(fmri_mat, step)

    # blank volume that will be used for ensemble average for fMRI volumes
    # prior to functional-anatomical coregistration
<<<<<<< .mine
    ave_fMRI_volume = np.zeros(fmri_desc['layers']*fmri_desc['rows']*fmri_desc['cols'], dtype=np.float64)
=======
    ave_fMRI_volume = np.zeros([fmri_desc['layers']*fmri_desc['rows']*fmri_desc['cols']],
        dtype=np.float64)
>>>>>>> .r4446

    count = 0
    number_volumes = len(fMRIdata)
    measures = np.zeros(number_volumes, dtype=metric_test)
    # load and perturb (rotation, translation) the fMRI volumes
    for i in fMRIdata:
        image = load_volume(fmri_desc, i)
        # random perturbation of angle, translation for each volume beyond the first
        if count == 0:
            fmri_series[count] = image
            count = count + 1
        else:
            x = np.random.random(6) - 0.5
            x = 10.0 * x
            fmri_series[count] = demo_rotate_fMRI_volume(image, fmri_desc, x)
            measures[count]['rotate'][0:6] = x[0:6]
            count = count + 1


    # load and register the fMRI volumes with volume_0 using normalized cross correlation metric
    imageF = fmri_series[0]
    if smooth_image:
        imageF = filter_image_3D(imageF, fmri_fwhm, ftype)
    for i in range(1, number_volumes):
        imageG = fmri_series[i]
        if smooth_image:
            imageG = filter_image_3D(imageG, fmri_fwhm, ftype)
        # the measure prior to alignment 
        measures[i]['cost'] = check_alignment(imageF, fmri_mat, imageG, fmri_mat, method='ncc',
                                              lite=histo_method, smhist=smooth_histo)
        x = register(imageF, fmri_mat, imageG, fmri_mat, lite=histo_method, method='ncc',
                       opt_method=optimizer_method, smhist=smooth_histo)
        measures[i]['align_rotate'][0:6] = x[0:6]
        measures[i]['align_cost'] = check_alignment(imageF, fmri_mat, imageG, fmri_mat,
                                                    method='ncc', lite=histo_method,
						    smhist=smooth_histo, alpha=x[0],
                                                    beta=x[1], gamma=x[2], Tx=x[3],
						    Ty=x[4], Tz=x[5])


    # align the volumes and average them for co-registration with the anatomical MRI 
    ave_fMRI_volume = fmri_series[0]['data'].astype(np.float64)
    for i in range(1, number_volumes):
        image = fmri_series[i]
        x[0:6] = measures[i]['align_rotate'][0:6]
        # overwrite the fMRI volume with the aligned volume
        fmri_series[i] = remap_image(image, x, resample='cubic')
        ave_fMRI_volume = ave_fMRI_volume + fmri_series[i]['data'].astype(np.float64)

    ave_fMRI_volume = (ave_fMRI_volume / float(number_volumes)).astype(np.uint8)
    ave_fMRI_volume = {'data' : ave_fMRI_volume, 'mat' : imageF['mat'], 
                       'dim' : imageF['dim'], 'fwhm' : imageF['fwhm']}
    # register (using normalized mutual information) with the anatomical MRI
    if smooth_image:
        imageF_anat = filter_image_3D(imageF_anat, anat_fwhm, ftype)

    x = register(imageF_anat, anat_mat, ave_fMRI_volume, fmri_mat, lite=histo_method,
                   method='nmi', opt_method=optimizer_method, smhist=smooth_histo)

    print 'functional-anatomical align parameters '
    print x
    for i in range(number_volumes):
        image = fmri_series[i]
        # overwrite the fMRI volume with the anatomical-aligned volume
        fmri_series[i] = remap_image(image, x, resample='cubic')

    return measures, imageF, fmri_series


def demo_fMRI_resample(imageF_anat, imageF_anat_mat, fmri_series):
    resampled_fmri_series = {}
    number_volumes = len(fmri_series)
    for i in range(number_volumes):
        resampled_fmri_series[i] = resize_image(fmri_series[i], imageF_anat_mat)

    return resampled_fmri_series


