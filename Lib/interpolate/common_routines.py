
from scipy_base import *
from scipy_base.fastumath import *

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        return asarray([a])
    elif type(a) is ArrayType and len(a)==1:
        # Takes care of mapping array(number) to array([number])
        return asarray([a[0]])
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



