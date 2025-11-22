import numpy as np
import scipy as scp

def spherical_jn_zeros(l, n_zeros):
	'''
 	 --------------------------------------------------------
	| Find the first n zeros of the spherical Bessel function |
 	 --------------------------------------------------------

	Parameters :
	--> l : int -- spherical Bessel function's order
	--> n_zeros : int -- number of roots to find

	Returns :
	--> zeros : ndarray -- array of roots from the spherical Bessel function

	Notes : 
	--> Roots are found by using Brent's method
	--> Each root is searched in an interval estimated by n : [(n+1/2)*pi;(1+n+1/2)*pi

	Examples :
	---------------

	>>> from scipy.special import spherical_jn_zeros
	>>> spherical_jn_zeros(0,10)
	array[ 3.14159265  6.28318531  9.42477796]

	---------------

	>>> from scipy.special import spherical_jn_zeros and spherical_jn
	>>> spherical_jn(0 ,spherical_jn_zeros(0,3))
	array[ 3.89804309e-17 -3.89804309e-17  3.89804309e-17]
	This test confirms the zeros, results are in an order of 10^(-17) so insignificant.

	---------------
	'''
	if not isinstance(l, int) or not isinstance(n_zeros, int) or l < 0 or n_zeros < 1:
		raise ValueError("`l` and `n_zeros` have to be natural integers")

	zeros = []
	for n in range(1, n_zeros + 1):
		x_min = (n + l/2) * np.pi
		x_max = (n + l/2 + 1) * np.pi
		root = scp.optimize.brentq(lambda x: scp.special.spherical_jn(l, x), x_min, x_max)
		zeros.append(root)
    
	return np.array(zeros)


