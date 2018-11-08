from . cimport c_zeros


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef scipy_zeros_parameters *myparams = <scipy_zeros_parameters *> params
    cdef void* args = myparams.args
    cdef callback_type f = myparams.function

    return f(x, args)  # callback_type takes a double and a struct


# cythonized way to call scalar bisect
cdef double bisect(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    myparams.root = c_zeros.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return myparams.root


# cythonized way to call scalar ridder
cdef double ridder(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    myparams.root = c_zeros.ridder(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return myparams.root


# cythonized way to call scalar brenth
cdef double brenth(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    myparams.root = c_zeros.brenth(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return myparams.root


# cythonized way to call scalar brentq
cdef double brentq(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    myparams.root = c_zeros.brentq(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return myparams.root


# Examples for testing
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
