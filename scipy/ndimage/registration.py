import math
import os
import numpy as N
import scipy.ndimage.register as R
import scipy.special  as SP
import scipy.ndimage  as NDI
import scipy.optimize as OPT
import time

# anatomical MRI to test with
# test registration on same image (optimal vector is (0,0,0,0,0,0)
inputname = 'ANAT1_V0001.img'
filename = os.path.join(os.path.split(__file__)[0], inputname)

def python_coreg(ftype=2, smimage=0, lite=1, smhist=0, method='mi', opt_method='powell'):
	# get_images is testing with 2 copies of anatomical MRI
	image1, image2, imdata = get_images(ftype, smimage)
    	start = time.time()
	parm_vector = multires_registration(image1, image2, imdata, lite, smhist, method, opt_method)
    	stop = time.time()
	print 'Total Optimizer Time is ', (stop-start)

	return parm_vector

def get_images(ftype, smimage):
	image1 = load_image()
	image2 = load_image()
	imdata = build_structs()
	image1['fwhm'] = build_fwhm(image1['mat'], imdata['step'])
	image2['fwhm'] = build_fwhm(image2['mat'], imdata['step'])
	if smimage: 
	    image_F_xyz1 = filter_image_3D(image1['data'], image1['fwhm'], ftype)
	    image1['data'] = image_F_xyz1
	    image_F_xyz2 = filter_image_3D(image2['data'], image2['fwhm'], ftype)
	    image2['data'] = image_F_xyz2

	return image1, image2, imdata

def multires_registration(image1, image2, imdata, lite, smhist, method, opt_method):
	ret_histo=0
	# make the step a scalar to can put in a multi-res loop
	loop = range(imdata['sample'].size)
    	x = imdata['parms']
	for i in loop:
	    step = imdata['sample'][i]
	    imdata['step'][:] = step
	    optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)
	    p_args = (optfunc_args,)
	    if opt_method=='powell':
		print 'POWELL multi-res registration step size ', step
		print 'vector ', x
    	        x = OPT.fmin_powell(optimize_function, x, args=p_args, callback=callback_powell) 
	    elif opt_method=='cg':
		print 'CG multi-res registration step size ', step
		print 'vector ', x
    	        x = OPT.fmin_cg(optimize_function, x, args=p_args, callback=callback_cg) 

	return x


def test_image_filter(image, imdata, ftype=2):
	image['fwhm'] = build_fwhm(image['mat'], imdata['step'])
	filt_image = filter_image_3D(image['data'], image['fwhm'], ftype)
	return filt_image

def test_optimizers(step=2, smooth=0, shist=0):
	opt_stats = {}
	print 'powell with stochastic resampling'
	x_0, p_time_0 = optimizer_powell(lite=0, smimage=smooth, smhist=shist, stepsize=step)
	opt_stats[0] = {'parms' : x_0, 'time' : p_time_0, 
			'label' : 'powell with stochastic resampling'}
	print 'powell without stochastic resampling'
	x_1, p_time_1 = optimizer_powell(lite=1, smimage=smooth, smhist=shist, stepsize=step)
	opt_stats[1] = {'parms' : x_1, 'time' : p_time_1, 
			'label' : 'powell without stochastic resampling'}
	print 'conjugate gradient with stochastic resampling'
	xcg_0, cg_time_0 = optimizer_cg(lite=0, smimage=smooth, smhist=shist, stepsize=step)
	opt_stats[2] = {'parms' : xcg_0, 'time' : cg_time_0, 
			'label' : 'conjugate gradient with stochastic resampling'}
	print 'conjugate gradient without stochastic resampling'
	xcg_1, cg_time_1 = optimizer_cg(lite=1, smimage=smooth, smhist=shist, stepsize=step)
	opt_stats[3] = {'parms' : xcg_1, 'time' : cg_time_1, 
			'label' : 'conjugate gradient without stochastic resampling'}
	return opt_stats

def callback_powell(x):
	print 'Parameter Vector from Powell: - '
	print x
	return

def callback_cg(x):
	print 'Parameter Vector from Conjugate Gradient: - '
	print x
	return

def optimizer_powell(lite=0, smhist=0, smimage=1, method='mi', ftype=2, stepsize=2):
	# test the Powell registration on the anatomical MRI volume
	image1 = load_image()
	image2 = load_image()
	imdata = build_structs(step=stepsize)
	image1['fwhm'] = build_fwhm(image1['mat'], imdata['step'])
	image2['fwhm'] = build_fwhm(image2['mat'], imdata['step'])
	M = build_rotate_matrix(imdata['parms'])
	if smimage: 
	    image_F_xyz1 = filter_image_3D(image1['data'], image1['fwhm'], ftype)
	    image1['data'] = image_F_xyz1
	    image_F_xyz2 = filter_image_3D(image2['data'], image2['fwhm'], ftype)
	    image2['data'] = image_F_xyz2

	ret_histo=0
	optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)
	p_args = (optfunc_args,)
    	start = time.time()
    	x = OPT.fmin_powell(optimize_function, imdata['parms'], args=p_args, callback=callback_powell) 
    	stop = time.time()
	return x, (stop-start)


def optimizer_cg(lite=0, smhist=0, smimage=1, method='mi', ftype=2, stepsize=2):
	# test the CG registration on the anatomical MRI volume
	image1 = load_image()
	image2 = load_image()
	imdata = build_structs(step=stepsize)
	image1['fwhm'] = build_fwhm(image1['mat'], imdata['step'])
	image2['fwhm'] = build_fwhm(image2['mat'], imdata['step'])
	M = build_rotate_matrix(imdata['parms'])
	if smimage: 
	    image_F_xyz1 = filter_image_3D(image1['data'], image1['fwhm'], ftype)
	    image1['data'] = image_F_xyz1
	    image_F_xyz2 = filter_image_3D(image2['data'], image2['fwhm'], ftype)
	    image2['data'] = image_F_xyz2

	ret_histo=0
	optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)
	p_args = (optfunc_args,)
    	start = time.time()
    	x = OPT.fmin_cg(optimize_function, imdata['parms'], args=p_args, callback=callback_cg) 
    	stop = time.time()
	return x, (stop-start)


def reg_single_pass(lite=0, smhist=0, smimage=0, method='mi', ftype=2, alpha=0.0, 
	            beta=0.0, gamma=0.0, Tx=0.0, Ty=0.0, Tz=0.0, ret_histo=0, stepsize=2):
	image1 = load_image()
	image2 = load_image()
	imdata = build_structs(step=stepsize)
	imdata['parms'][0] = alpha
	imdata['parms'][1] = beta
	imdata['parms'][2] = gamma
	imdata['parms'][3] = Tx
	imdata['parms'][4] = Ty
	imdata['parms'][5] = Tz
	image1['fwhm'] = build_fwhm(image1['mat'], imdata['step'])
	image2['fwhm'] = build_fwhm(image2['mat'], imdata['step'])
	print 'image1[fwhm] ', image1['fwhm'] 
	print 'image2[fwhm] ', image2['fwhm'] 
	M = build_rotate_matrix(imdata['parms'])
	if smimage: 
	    image_F_xyz1 = filter_image_3D(image1['data'], image1['fwhm'], ftype)
	    image1['data'] = image_F_xyz1
	    image_F_xyz2 = filter_image_3D(image2['data'], image2['fwhm'], ftype)
	    image2['data'] = image_F_xyz2

	optfunc_args = (image1, image2, imdata['step'], imdata['fwhm'], lite, smhist, method, ret_histo)

	if ret_histo:
	    cost, joint_histogram = optimize_function(imdata['parms'], optfunc_args)
	    return cost, joint_histogram 
    	else:
	    cost = optimize_function(imdata['parms'], optfunc_args)
	    return cost


def smooth_kernel(fwhm, x, ktype=1):
	eps = 0.00001
	s   = N.square((fwhm/math.sqrt(8.0*math.log(2.0)))) + eps
	if ktype==1:
	    # from SPM: Gauss kernel convolved with 1st degree B spline
	    w1 = 0.5 * math.sqrt(2.0/s)
	    w2 = -0.5 / s
	    w3 = math.sqrt((s*math.pi) /2.0)
	    kernel = 0.5*(SP.erf(w1*(x+1))*(x+1)       + SP.erf(w1*(x-1))*(x-1)    - 2.0*SP.erf(w1*(x))*(x) + 
	 	          w3*(N.exp(w2*N.square(x+1))) + N.exp(w2*(N.square(x-1))) - 2.0*N.exp(w2*N.square(x)))
	    kernel[kernel<0] = 0
	    kernel = kernel / kernel.sum()  
	else:
	    # Gauss kernel 
	    kernel = (1.0/math.sqrt(2.0*math.pi*s)) * N.exp(-N.square(x)/(2.0*s)) 
	    kernel = kernel / kernel.sum()  

	return kernel

def filter_image_3D(imageRaw, fwhm, ftype=2):
	p = N.ceil(2*fwhm[0]).astype(int)
	x = N.array(range(-p, p+1))
	kernel_x = smooth_kernel(fwhm[0], x, ktype=ftype)
	p = N.ceil(2*fwhm[1]).astype(int)
	x = N.array(range(-p, p+1))
	kernel_y = smooth_kernel(fwhm[1], x, ktype=ftype)
	p = N.ceil(2*fwhm[2]).astype(int)
	x = N.array(range(-p, p+1))
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
	view_3x3 = N.square(M[0:3, 0:3])
	vxg = N.sqrt(view_3x3.sum(axis=0))
	# assumes that sampling is the same for xyz
 	size = N.array([1,1,1])*S[0]
	x = N.square(size) - N.square(vxg)
	# clip
	x[x<0] = 0
	fwhm = N.sqrt(x) / vxg
	# pathology when stepsize = 1 for MAT equal to the identity matrix
	fwhm[fwhm==0] = 1
	# return the 3D Gaussian kernel width (xyz)
	return fwhm 

def load_image(imagename=filename, rows=256, cols=256, layers=90):
    	ImageVolume = N.fromfile(imagename, dtype=N.uint16).reshape(layers, rows, cols);
	# clip to 8 bits. this will be rescale to 8 bits for fMRI
    	ImageVolume[ImageVolume>255] = 255
	# voxel to pixel is identity for this simulation using anatomical MRI volume
	# 4x4 matrix
	M = N.eye(4, dtype=N.float64);
	# dimensions 
	D = N.zeros(3, dtype=N.int32);
	# Gaussian kernel - fill in with build_fwhm() 
	F = N.zeros(3, dtype=N.float64);
	D[0] = rows
	D[1] = cols
	D[2] = layers
	# make sure the data type is uchar
    	ImageVolume = ImageVolume.astype(N.uint8)
	image = {'data' : ImageVolume, 'mat' : M, 'dim' : D, 'fwhm' : F}
    	return image


def optimize_function(x, optfunc_args):
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
	# image_F is base image
	# image_G is the rotated image
	# rot_matrix is the 4x4 constructed (current angles and translates) transform matrix
	# sample_vector is the subsample vector for x-y-z

	F_inv = N.linalg.inv(image_F['mat'])
	composite = N.dot(F_inv, rot_matrix)
	composite = N.dot(composite, image_G['mat'])

	# allocate memory from Python as memory leaks when created in C-ext
	joint_histogram = N.zeros([256, 256], dtype=N.float64);

	if do_lite: 
	    R.register_histogram_lite(image_F['data'], image_G['data'], composite, sample_vector, joint_histogram)
	else:
	    R.register_histogram(image_F['data'], image_G['data'], composite, sample_vector, joint_histogram)

	if smooth: 
	    p = N.ceil(2*fwhm[0]).astype(int)
	    x = N.array(range(-p, p+1))
	    kernel1 = smooth_kernel(fwhm[0], x)
	    p = N.ceil(2*fwhm[1]).astype(int)
	    x = N.array(range(-p, p+1))
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
	    marginal_outer = N.outer(marginal_col, marginal_row)
	    H = joint_histogram * N.log(joint_histogram / marginal_outer)  
	    mutual_information = H.sum()
	    cost = -mutual_information

	elif method == 'ecc':
	    # entropy correlation coefficient 
	    marginal_outer = N.outer(marginal_col, marginal_row)
	    H = joint_histogram * N.log(joint_histogram / marginal_outer)  
	    mutual_information = H.sum()
	    row_entropy = marginal_row * N.log(marginal_row)
	    col_entropy = marginal_col * N.log(marginal_col)
	    ecc  = -2.0*mutual_information/(row_entropy.sum() + col_entropy.sum())
	    cost = -ecc

	elif method == 'nmi':
	    # normalized mutual information
	    row_entropy = marginal_row * N.log(marginal_row)
	    col_entropy = marginal_col * N.log(marginal_col)
	    H = joint_histogram * N.log(joint_histogram)  
	    nmi = (row_entropy.sum() + col_entropy.sum()) / (H.sum())
	    cost = -nmi

	elif method == 'ncc':
	    # cross correlation from the joint histogram 
	    r, c = joint_histogram.shape
	    i = N.array(range(1,c+1))
	    j = N.array(range(1,r+1))
	    m1 = (marginal_row * i).sum()
	    m2 = (marginal_col * j).sum()
	    sig1 = N.sqrt((marginal_row*(N.square(i-m1))).sum())
	    sig2 = N.sqrt((marginal_col*(N.square(j-m2))).sum())
	    [a, b] = N.mgrid[1:c+1, 1:r+1]
	    a = a - m1
	    b = b - m2
	    # element multiplies in the joint histogram and grids
	    H = ((joint_histogram * a) * b).sum()
	    ncc = H / (N.dot(sig1, sig2)) 
	    cost = -ncc

	if ret_histo:
	    return cost, joint_histogram 
    	else:
	    return cost


def build_structs(step=2):
	# build image data structures here
	P = N.zeros(6, dtype=N.float64);
	T = N.zeros(6, dtype=N.float64);
	F = N.zeros(2, dtype=N.int32);
	S = N.ones(3,  dtype=N.int32);
	sample = N.zeros(2, dtype=N.int32);
	S[0] = step
	S[1] = step
	S[2] = step
	# histogram smoothing
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
	# P[0] = alpha <=> pitch
	# P[1] = beta  <=> roll
	# P[2] = gamma <=> yaw
	# P[3] = Tx
	# P[4] = Ty
	# P[5] = Tz
	img_data = {'parms' : P, 'step' : S, 'fwhm' : F, 'tol' : T, 'sample' : sample}
	return img_data


def build_rotate_matrix(img_data_parms):
	R1 = N.zeros([4,4], dtype=N.float64);
	R2 = N.zeros([4,4], dtype=N.float64);
	R3 = N.zeros([4,4], dtype=N.float64);
	T  = N.eye(4, dtype=N.float64);

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

	rot_matrix = N.dot(T, R1);
	rot_matrix = N.dot(rot_matrix, R2);
	rot_matrix = N.dot(rot_matrix, R3);

	return rot_matrix






