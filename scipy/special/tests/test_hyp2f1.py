from __future__ import division, print_function, absolute_import

import itertools
import numpy as np

from scipy.special._testutils import FuncData, MissingModule, check_version
from scipy.special import hyp2f1

try:
    import mpmath
except ImportError:
    mpmath = MissingModule('mpmath')

# Fix the seed number
np.random.seed(1108257)

def _get_rand(minv, maxv, only_int=False, exclude_int=False):
    '''Gets a random value in the range (minv, maxv).
    
    When only_int=True, this function only returns integers in (minv, maxv).
    When exclude_int=True, this function only returns *non*-integers in (minv, maxv).
    '''
    if only_int:
        randv = np.random.randint(minv,maxv+1)
    else:
        keep_going = True
        count = 0
        while keep_going and count < 100:
            randv = minv + (maxv - minv)*np.random.rand()
            # Get a new randv if we are excluding integers and this value is too close to an int
            keep_going = exclude_int and np.abs(np.round(randv)-randv) < 1.e-2
            count += 1

    return randv

def _get_rnd_tuple(*args):
    '''
    Returns a tuple containing randomly generated values.

    The returned tuple will contain as many random values as there are arguments
    to this function. Each argument to _get_rnd_tuple should be a 3-tuple of the form:

    (flag, minv, maxv)

    flag: i -> this value will be a random integer
            f -> this value will be a floating point number BUT not an integer
            a -> this value will be any random number between minv, maxv
                    (integers included)
    '''
    curr_params = []

    for p in args:
        flag = p[0]
        flag = flag.capitalize()

        minv = p[1]
        maxv = p[2]

        # I -> pick a random integer
        if flag == 'I':
            randv = _get_rand(minv, maxv, only_int=True)
        # F -> pick a floating point number (exclude integers)
        elif flag == 'F':
            randv = _get_rand(minv, maxv, exclude_int=True)
        else:
            randv = _get_rand(minv, maxv)

        curr_params.append(randv)

    return tuple(curr_params)

def _add_params(plist,*args):
    rnd_tuple = _get_rnd_tuple(*args)
    plist.append(rnd_tuple)
    return

# Build a list of a, b, c values that covers each code branch
# in specfun.f
def _build_abc_list():
    small = 1.e-2
    plist = []

    # +a, +b, +c test case
    _add_params(plist,('f',small,3.0),('f',small,3.0),('f',small,3.0))

    # -a, +b, +c test case
    _add_params(plist,('f',-3.0,-small),('f',small,3.0),('f',small,3.0))

    # -a, -b, +c test case
    _add_params(plist,('f',-3.0,-small),('f',-3.0,-small),('f',small,3.0))

    # -a, -b, -c test case
    _add_params(plist,('f',-3.0,-small),('f',-3.0,-small),('f',-3.0,-small))

    # Re(c-a-b)>0 test case
    _add_params(plist,('f',small,2.0),('f',small,2.0),('f',4.0+small,6.0))

    # Re(c-a-b)<-1 test case
    _add_params(plist,('f',2.0,4.0),('f',2.0,4.0),('f',small,3.0-small))

    # c-a-b=m>0
    a, b, m = _get_rnd_tuple(('f',small,2.0),('f',small,2.0),('i',1,6))
    c = a + b + m
    plist.append((a,b,c))

    # c-a-b=0
    a, b = _get_rnd_tuple(('f',small,2.0),('f',small,2.0))
    c = a + b
    plist.append((a,b,c))

    # c-a-b=-m<0
    a, b, m = _get_rnd_tuple(('f',3.0,5.0),('f',3.0,5.0),('i',1,6))
    c = a + b - m
    plist.append((a,b,c))

    # b-a=m>0
    a, m, c = _get_rnd_tuple(('f',small,2.0),('i',1,5),('f',small,2.0))
    b = a + m
    plist.append((a,b,c))

    # b-a=0
    a, c = _get_rnd_tuple(('f',small,3.0),('f',small,3.0))
    b = a
    plist.append((a,b,c))

    # b-a=-m<0
    a, m, c = _get_rnd_tuple(('f',4.0,6.0),('i',1,4),('f',small,2.0))
    b = a - m
    plist.append((a,b,c))

    # c-a=m>0
    a, b, m = _get_rnd_tuple(('f',small,3.0),('f',small,3.0),('i',1,5))
    c = a + m
    plist.append((a,b,c))

    # c-a=0
    a, b = _get_rnd_tuple(('f',small,3.0),('f',small,3.0))
    c = a
    plist.append((a,b,c))
    
    # c-a=-m<0
    a, b, m = _get_rnd_tuple(('f',4.0,6.0),('f',small,3.0),('i',1,4))
    c = a - m
    plist.append((a,b,c))
    
    # a=-m
    m, b, c = _get_rnd_tuple(('i',1,10),('f',small,3.0),('f',small,3.0))
    a = -m
    plist.append((a,b,c))

    # c=-m, a=-n, |m|>|n|
    #b, m, n = _get_rnd_tuple(('f',small,3.0),('i',5,8),('i',2,4))
    #a = -n
    #c = -m
    #plist.append((a,b,c))
    
    return plist

# Build the test case data that FuncData needs
def _build_test_cases():
    abc_list = _build_abc_list()

    # Check each region of the complex plane that
    # uses a different method to calculate hyp2f1
    rho = [0.5, 0.95, 1.0, 1.05, 10.0]
    phi = np.pi*np.array([0.0, 0.1, -0.1, 0.375, -0.375, 0.5, -0.5, -1.0, -0.75, 0.75],dtype=np.float)

    # The total number of test cases
    N = 2*len(abc_list)*len(rho)*len(phi)
    dataset = np.zeros((N,5),dtype=np.complex)

    count = 0
    for (a, b, c), r, theta in itertools.product(abc_list,rho,phi):
        z = r*np.exp(1.j*theta)

        # For most values of a,b,c mpmath returns the x - 0j branch
        # of hyp2f1 on the branch cut x=(1,inf) whereas scipy's
        # hyp2f1 calculates the x + 0j branch. Thus, to generate the right
        # comparison values on the branch cut, we evaluate mpmath.hyp2f1
        # at x + 1e-15*j.
        #
        # The exception to this occurs when c-a=-m in which case both
        # mpmath and scipy calculate the x + 0j branch on the branch
        # cut. When this happens mpmath.hyp2f1 will be evaluated
        # at the normal z point.
        on_branch_cut = (np.real(z) > 1.0) and (np.abs(np.imag(z)) < 1.e-15)

        i_ca = np.round(c-a)
        i_cb = np.round(c-b)

        c_a_int = (np.abs(c-a-i_ca) < 1.e-15) and (i_ca <= 0)
        c_b_int = (np.abs(c-b-i_cb) < 1.e-15) and (i_cb <= 0)

        # Make sure imaginary part is *exactly* zero
        if on_branch_cut:
            z = np.real(z) + 0.0j
        
        if on_branch_cut and not (c_a_int or c_b_int):
            z_mpmath = np.real(z) + 1.e-15j
        else:
            z_mpmath = z

        right = np.complex(mpmath.hyp2f1(a, b, c, z_mpmath))
        dataset[count,:] = a, b, c, z, right
        count += 1

        # Make sure that scipy's hyp2f1 returns same result
        # when we swap a, b
        dataset[count,:] = b, a, c, z, right
        count += 1
    
    return dataset

def _hyp2f1_wrap(a,b,c,z):
    a = np.real(a)
    b = np.real(b)
    c = np.real(c)
    return hyp2f1(a,b,c,z)

##########################################
#  Testing hyp2f1 throughout the complex plane
##########################################
@check_version(mpmath, '1.0.0')
def test_hyp2f1():
    dataset = _build_test_cases()

    FuncData(_hyp2f1_wrap, dataset, (0,1,2,3j), 4, rtol=1.e-5).check()
    return
