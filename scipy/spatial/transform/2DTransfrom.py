#2D matrices used for linear transformation in linear algebra 

import numpy as np

def scaling(scaling_factor):
	return(scaling_factor*(np.identity(2))) #This returns a scaling matrix

def rotation(angle):
	c = np.cos(angle)
	s = np.sin(angle)
	arr = np.empty(c.shape + (2, 2), dtype=c.dtype)
	arr[..., 0, 0] = c
	arr[..., 0, 1] = -s
	arr[..., 1, 0] = s
	arr[..., 1, 1] = c

	return arr #This returns a rotation matrix

	
def projection(angle):
	c = np.cos(angle)
	s = np.sin(angle)
	arr = np.empty(c.shape + (2, 2), dtype=c.dtype)
	arr[..., 0, 0] = c*c
	arr[..., 0, 1] = c*s
	arr[..., 1, 0] = c*s
	arr[..., 1, 1] = s*s

	return arr #This returns a rotation matrix


def reflection(angle):
	arr = np.empty(c.shape + (2, 2), dtype=c.dtype)
	c = np.cos(angle)
	s = np.sin(angle)
	arr[..., 0, 0] = (2*c) -1
	arr[..., 0, 1] = 2*c*s
	arr[..., 1, 0] = 2*c*s
	arr[..., 1, 1] = (2*s)-1

	return arr #This returns a reflection matrix
