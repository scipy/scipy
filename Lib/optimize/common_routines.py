from scipy_base import *
from scipy_base.fastumath import *

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        a = asarray([a])
    if len(a) == 0:
        a = asarray([a])
    return a

def check_func(thefunc, x0, args, numinputs, output_shape=None):
    args = (x0[:numinputs],) + args
    res = myasarray(apply(thefunc,args))
    if (output_shape != None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            raise TypeError, "There is a mismatch between the input and output shape of %s." % thefunc.func_name
    return shape(res)



