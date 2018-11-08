"""
Use a structure to hold callback args so there's no Python and functions can
free the global interpreter lock.
"""

from ...cython_optimize cimport (
    brentq, brenth, ridder, bisect, scipy_zeros_parameters)

ARGS = (0.0, 0.0, 1.0)
A0 = tuple(-2.0 - x/10.0 for x in range(10))
XLO, XHI = 0.0, 2.0
XTOL, RTOL, MITR = 0.001, 0.001, 10


# extra parameters
ctypedef struct extra_params:
    double[4] a


# full output structure
ctypedef struct full_output_struct:
    int funcalls
    int iterations
    int error_num
    double root


# callback function
cdef double f_example(double x, void *args):
    cdef extra_params *myargs = <extra_params *> args
    cdef double[4] a
    # unpack structure
    a = myargs.a
    return a[3]*x**3 + a[2]*x**2 + a[1]*x + a[0]


# bisect example
cdef double bisect_example(tuple args):
    cdef extra_params myargs
    myargs.a = args
    return bisect(
        f_example, XLO, XHI, &myargs, XTOL, RTOL, MITR, NULL)


# ridder example
cdef double ridder_example(tuple args):
    cdef extra_params myargs
    myargs.a = args
    return ridder(
        f_example, XLO, XHI, &myargs, XTOL, RTOL, MITR, NULL)


# brenth example
cdef double brenth_example(tuple args):
    cdef extra_params myargs
    myargs.a = args
    return brenth(
        f_example, XLO, XHI, &myargs, XTOL, RTOL, MITR, NULL)


# brentq example
cdef double brentq_example(tuple args):
    cdef extra_params myargs
    myargs.a = args
    return brentq(
        f_example, XLO, XHI, &myargs, XTOL, RTOL, MITR, NULL)


# python function
def loop_example(method, a0=A0, args=ARGS):
    """example of zeros functions in a loop"""
    method = method.lower()
    args = [(a0_,) + args for a0_ in a0]
    if method == 'bisect':
        return map(bisect_example, args)
    elif method == 'ridder':
        return map(ridder_example, args)
    elif method == 'brenth':
        return map(brenth_example, args)
    elif method == 'brentq':
        return map(brentq_example, args)
    else:
        raise ValueError('unrecognized method')


# brentq example with full ouptut
cdef full_output_struct brentq_full_output_example(tuple args):
    cdef full_output_struct my_full_output
    cdef scipy_zeros_parameters full_output
    cdef extra_params myargs
    myargs.a = args
    my_full_output.root = brentq(
        f_example, XLO, XHI, &myargs, XTOL, RTOL, MITR, &full_output)
    my_full_output.funcalls = full_output.funcalls
    my_full_output.iterations = full_output.iterations
    my_full_output.error_num = full_output.error_num
    return my_full_output


# python function
def full_output_example():
    return brentq_full_output_example((A0[0],) + ARGS)
