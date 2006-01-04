## Automatically adapted for scipy Oct 21, 2005 by 


from numpy import *
from numpy.umath import *

def myasarray(a):
    if isscalar(a):
        return asarray([a])
    else:
        return asarray(a)

def check_func(thefunc, x0, args, numinputs, output_shape=None):
    args = (x0[:numinputs],) + args
    res = myasarray(apply(thefunc,args))
    if (output_shape != None) and (res.shape != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return res.shape
            raise TypeError, "There is a mismatch between the input and output shape of %s." % thefunc.func_name
    return res.shape



