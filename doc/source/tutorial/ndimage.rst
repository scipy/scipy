Multi-dimensional image processing (`scipy.ndimage`)
====================================================

.. moduleauthor:: Peter Verveer <verveer@users.sourceforge.net>

.. currentmodule:: scipy.ndimage

.. _ndimage-introduction:

Introduction
------------

Image processing and analysis are generally seen as operations on
two-dimensional arrays of values. There are however a number of
fields where images of higher dimensionality must be analyzed. Good
examples of these are medical imaging and biological imaging.
:mod:`numpy` is suited very well for this type of applications due
its inherent multi-dimensional nature. The :mod:`scipy.ndimage`
packages provides a number of general image processing and analysis
functions that are designed to operate with arrays of arbitrary
dimensionality. The packages currently includes functions for
linear and non-linear filtering, binary morphology, B-spline
interpolation, and object measurements.

.. _ndimage-properties-shared-by-all-functions:

Properties shared by all functions
----------------------------------

All functions share some common properties. Notably, all functions
allow the specification of an output array with the *output*
argument. With this argument you can specify an array that will be
changed in-place with the result with the operation. In this case
the result is not returned. Usually, using the *output* argument is
more efficient, since an existing array is used to store the
result.

The type of arrays returned is dependent on the type of operation,
but it is in most cases equal to the type of the input. If,
however, the *output* argument is used, the type of the result is
equal to the type of the specified output argument. If no output
argument is given, it is still possible to specify what the result
of the output should be. This is done by simply assigning the
desired `numpy` type object to the output argument. For example:

::

    >>> correlate(np.arange(10), [1, 2.5])
    array([ 0,  2,  6,  9, 13, 16, 20, 23, 27, 30])
    >>> correlate(np.arange(10), [1, 2.5], output=np.float64)
    array([  0. ,   2.5,   6. ,   9.5,  13. ,  16.5,  20. ,  23.5,  27. ,  30.5])

.. _ndimage-filter-functions:

Filter functions
----------------

.. currentmodule:: scipy.ndimage.filters

The functions described in this section all perform some type of spatial
filtering of the the input array: the elements in the output are some function
of the values in the neighborhood of the corresponding input element. We refer
to this neighborhood of elements as the filter kernel, which is often
rectangular in shape but may also have an arbitrary footprint. Many
of the functions described below allow you to define the footprint
of the kernel, by passing a mask through the *footprint* parameter.
For example a cross shaped kernel can be defined as follows:

::

    >>> footprint = array([[0,1,0],[1,1,1],[0,1,0]])
    >>> footprint
    array([[0, 1, 0],
           [1, 1, 1],
           [0, 1, 0]])

Usually the origin of the kernel is at the center calculated by
dividing the dimensions of the kernel shape by two. For instance,
the origin of a one-dimensional kernel of length three is at the
second element. Take for example the correlation of a
one-dimensional array with a filter of length 3 consisting of
ones:

::

    >>> a = [0, 0, 0, 1, 0, 0, 0]
    >>> correlate1d(a, [1, 1, 1])
    array([0, 0, 1, 1, 1, 0, 0])

Sometimes it is convenient to choose a different origin for the
kernel. For this reason most functions support the *origin*
parameter which gives the origin of the filter relative to its
center. For example:

::

    >>> a = [0, 0, 0, 1, 0, 0, 0]
    >>> correlate1d(a, [1, 1, 1], origin = -1)
    array([0 1 1 1 0 0 0])

The effect is a shift of the result towards the left. This feature
will not be needed very often, but it may be useful especially for
filters that have an even size. A good example is the calculation
of backward and forward differences:

::

    >>> a = [0, 0, 1, 1, 1, 0, 0]
    >>> correlate1d(a, [-1, 1])               # backward difference
    array([ 0  0  1  0  0 -1  0])
    >>> correlate1d(a, [-1, 1], origin = -1)  # forward difference
    array([ 0  1  0  0 -1  0  0])

We could also have calculated the forward difference as follows:

::

    >>> correlate1d(a, [0, -1, 1])
    array([ 0  1  0  0 -1  0  0])

However, using the origin parameter instead of a larger kernel is
more efficient. For multi-dimensional kernels *origin* can be a
number, in which case the origin is assumed to be equal along all
axes, or a sequence giving the origin along each axis.

Since the output elements are a function of elements in the
neighborhood of the input elements, the borders of the array need
to be dealt with appropriately by providing the values outside the
borders. This is done by assuming that the arrays are extended
beyond their boundaries according certain boundary conditions. In
the functions described below, the boundary conditions can be
selected using the *mode* parameter which must be a string with the
name of the boundary condition. Following boundary conditions are
currently supported:

 ==========   ====================================   ====================
 "nearest"    Use the value at the boundary          [1 2 3]->[1 1 2 3 3]
 "wrap"       Periodically replicate the array       [1 2 3]->[3 1 2 3 1]
 "reflect"    Reflect the array at the boundary      [1 2 3]->[1 1 2 3 3]
 "constant"   Use a constant value, default is 0.0   [1 2 3]->[0 1 2 3 0]
 ==========   ====================================   ====================

The "constant" mode is special since it needs an additional
parameter to specify the constant value that should be used.

.. note:: The easiest way to implement such boundary conditions would be to copy the data to a larger array and extend the data at the borders according to the boundary conditions. For large arrays and large filter kernels, this would be very memory consuming, and the functions described below therefore use a different approach that does not require allocating large temporary buffers.

Correlation and convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The :func:`correlate1d` function calculates a one-dimensional correlation
    along the given axis. The lines of the array along the given axis
    are correlated with the given *weights*. The *weights* parameter
    must be a one-dimensional sequences of numbers.


    The function :func:`correlate` implements multi-dimensional correlation
    of the input array with a given kernel.


    The :func:`convolve1d` function calculates a one-dimensional convolution
    along the given axis. The lines of the array along the given axis
    are convoluted with the given *weights*. The *weights* parameter
    must be a one-dimensional sequences of numbers.

    .. note:: A convolution is essentially a correlation after mirroring the kernel. As a result, the *origin* parameter behaves differently than in the case of a correlation: the result is shifted in the opposite directions.

    The function :func:`convolve` implements multi-dimensional convolution of
    the input array with a given kernel.

    .. note:: A convolution is essentially a correlation after mirroring the kernel. As a result, the *origin* parameter behaves differently than in the case of a correlation: the results is shifted in the opposite direction.

.. _ndimage-filter-functions-smoothing:

Smoothing filters
^^^^^^^^^^^^^^^^^


    The :func:`gaussian_filter1d` function implements a one-dimensional
    Gaussian filter. The standard-deviation of the Gaussian filter is
    passed through the parameter *sigma*. Setting *order* = 0 corresponds
    to convolution with a Gaussian kernel. An order of 1, 2, or 3
    corresponds to convolution with the first, second or third
    derivatives of a Gaussian. Higher order derivatives are not
    implemented.


    The :func:`gaussian_filter` function implements a multi-dimensional
    Gaussian filter. The standard-deviations of the Gaussian filter
    along each axis are passed through the parameter *sigma* as a
    sequence or numbers. If *sigma* is not a sequence but a single
    number, the standard deviation of the filter is equal along all
    directions. The order of the filter can be specified separately for
    each axis. An order of 0 corresponds to convolution with a Gaussian
    kernel. An order of 1, 2, or 3 corresponds to convolution with the
    first, second or third derivatives of a Gaussian. Higher order
    derivatives are not implemented. The *order* parameter must be a
    number, to specify the same order for all axes, or a sequence of
    numbers to specify a different order for each axis.

    .. note:: The multi-dimensional filter is implemented as a sequence of one-dimensional Gaussian filters. The intermediate arrays are stored in  the same data type as the output.  Therefore, for output types with a lower precision, the results may be imprecise because intermediate results may be stored with insufficient precision. This can be prevented by specifying a more precise output type.


    The :func:`uniform_filter1d` function calculates a one-dimensional
    uniform filter of the given *size* along the given axis.


    The :func:`uniform_filter` implements a multi-dimensional uniform
    filter. The sizes of the uniform filter are given for each axis as
    a sequence of integers by the *size* parameter. If *size* is not a
    sequence, but a single number, the sizes along all axis are assumed
    to be equal.

    .. note:: The multi-dimensional filter is implemented as a sequence of one-dimensional uniform filters. The intermediate arrays are stored in the same data type as the output. Therefore, for output types with a lower precision, the results may be imprecise because intermediate results may be stored with insufficient precision. This can be prevented by specifying a more precise output type.


Filters based on order statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The :func:`minimum_filter1d` function calculates a one-dimensional
    minimum filter of given *size* along the given axis.


    The :func:`maximum_filter1d` function calculates a one-dimensional
    maximum filter of given *size* along the given axis.


    The :func:`minimum_filter` function calculates a multi-dimensional
    minimum filter. Either the sizes of a rectangular kernel or the
    footprint of the kernel must be provided. The *size* parameter, if
    provided, must be a sequence of sizes or a single number in which
    case the size of the filter is assumed to be equal along each axis.
    The *footprint*, if provided, must be an array that defines the
    shape of the kernel by its non-zero elements.


    The :func:`maximum_filter` function calculates a multi-dimensional
    maximum filter. Either the sizes of a rectangular kernel or the
    footprint of the kernel must be provided. The *size* parameter, if
    provided, must be a sequence of sizes or a single number in which
    case the size of the filter is assumed to be equal along each axis.
    The *footprint*, if provided, must be an array that defines the
    shape of the kernel by its non-zero elements.


    The :func:`rank_filter` function calculates a multi-dimensional rank
    filter. The *rank* may be less then zero, i.e., *rank* = -1 indicates
    the largest element. Either the sizes of a rectangular kernel or
    the footprint of the kernel must be provided. The *size* parameter,
    if provided, must be a sequence of sizes or a single number in
    which case the size of the filter is assumed to be equal along each
    axis. The *footprint*, if provided, must be an array that defines
    the shape of the kernel by its non-zero elements.


    The :func:`percentile_filter` function calculates a multi-dimensional
    percentile filter. The *percentile* may be less then zero, i.e.,
    *percentile* = -20 equals *percentile* = 80. Either the sizes of a
    rectangular kernel or the footprint of the kernel must be provided.
    The *size* parameter, if provided, must be a sequence of sizes or a
    single number in which case the size of the filter is assumed to be
    equal along each axis. The *footprint*, if provided, must be an
    array that defines the shape of the kernel by its non-zero
    elements.


    The :func:`median_filter` function calculates a multi-dimensional median
    filter. Either the sizes of a rectangular kernel or the footprint
    of the kernel must be provided. The *size* parameter, if provided,
    must be a sequence of sizes or a single number in which case the
    size of the filter is assumed to be equal along each axis. The
    *footprint* if provided, must be an array that defines the shape of
    the kernel by its non-zero elements.


Derivatives
^^^^^^^^^^^

Derivative filters can be constructed in several ways. The function
:func:`gaussian_filter1d` described in  
:ref:`ndimage-filter-functions-smoothing` can be used to calculate
derivatives along a given axis using the *order* parameter. Other
derivative filters are the Prewitt and Sobel filters:

    The :func:`prewitt` function calculates a derivative along the given
    axis.


    The :func:`sobel` function calculates a derivative along the given
    axis.


The Laplace filter is calculated by the sum of the second
derivatives along all axes. Thus, different Laplace filters can be
constructed using different second derivative functions. Therefore
we provide a general function that takes a function argument to
calculate the second derivative along a given direction and to
construct the Laplace filter:

    The function :func:`generic_laplace` calculates a laplace filter using
    the function passed through :func:`derivative2` to calculate second
    derivatives. The function :func:`derivative2` should have the following
    signature::

         derivative2(input, axis, output, mode, cval, *extra_arguments, **extra_keywords)
         
    It should calculate the second derivative along the dimension
    *axis*. If *output* is not None it should use that for the output
    and return None, otherwise it should return the result. *mode*,
    *cval* have the usual meaning.

    The *extra_arguments* and *extra_keywords* arguments can be used
    to pass a tuple of extra arguments and a dictionary of named
    arguments that are passed to :func:`derivative2` at each call.

    For example::

        >>> def d2(input, axis, output, mode, cval):
        ...     return correlate1d(input, [1, -2, 1], axis, output, mode, cval, 0)
        ... 
        >>> a = zeros((5, 5))
        >>> a[2, 2] = 1
        >>> generic_laplace(a, d2)
        array([[ 0.,  0.,  0.,  0.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  1., -4.,  1.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  0.,  0.,  0.,  0.]])

    To demonstrate the use of the *extra_arguments* argument we could
    do::

        >>> def d2(input, axis, output, mode, cval, weights):
        ...     return correlate1d(input, weights, axis, output, mode, cval, 0,)
        ... 
        >>> a = zeros((5, 5))
        >>> a[2, 2] = 1
        >>> generic_laplace(a, d2, extra_arguments = ([1, -2, 1],))
        array([[ 0.,  0.,  0.,  0.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  1., -4.,  1.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  0.,  0.,  0.,  0.]])

    or::

        >>> generic_laplace(a, d2, extra_keywords = {'weights': [1, -2, 1]})
        array([[ 0.,  0.,  0.,  0.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  1., -4.,  1.,  0.],
               [ 0.,  0.,  1.,  0.,  0.],
               [ 0.,  0.,  0.,  0.,  0.]])


The following two functions are implemented using
:func:`generic_laplace` by providing appropriate functions for the
second derivative function:

    The function :func:`laplace` calculates the Laplace using discrete
    differentiation for the second derivative (i.e. convolution with
    :obj:`[1, -2, 1]`).


    The function :func:`gaussian_laplace` calculates the Laplace using
    :func:`gaussian_filter` to calculate the second derivatives. The
    standard-deviations of the Gaussian filter along each axis are
    passed through the parameter *sigma* as a sequence or numbers. If
    *sigma* is not a sequence but a single number, the standard
    deviation of the filter is equal along all directions.


The gradient magnitude is defined as the square root of the sum of
the squares of the gradients in all directions. Similar to the
generic Laplace function there is a :func:`generic_gradient_magnitude`
function that calculated the gradient magnitude of an array:

    The function :func:`generic_gradient_magnitude` calculates a gradient
    magnitude using the function passed through :func:`derivative` to
    calculate first derivatives. The function :func:`derivative` should have
    the following signature::

        derivative(input, axis, output, mode, cval, *extra_arguments, **extra_keywords)
        
    
    It should calculate the derivative along the dimension *axis*. If
    *output* is not None it should use that for the output and return
    None, otherwise it should return the result. *mode*, *cval* have
    the usual meaning.

    The *extra_arguments* and *extra_keywords* arguments can be used
    to pass a tuple of extra arguments and a dictionary of named
    arguments that are passed to *derivative* at each call.

    For example, the :func:`sobel` function fits the required signature::

        >>> a = zeros((5, 5))
        >>> a[2, 2] = 1
        >>> generic_gradient_magnitude(a, sobel)
        array([[ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
               [ 0.        ,  1.41421356,  2.        ,  1.41421356,  0.        ],
               [ 0.        ,  2.        ,  0.        ,  2.        ,  0.        ],
               [ 0.        ,  1.41421356,  2.        ,  1.41421356,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])

    See the documentation of :func:`generic_laplace` for examples of using
    the *extra_arguments* and *extra_keywords* arguments.


The :func:`sobel` and :func:`prewitt` functions fit the required signature and
can therefore directly be used with :func:`generic_gradient_magnitude`.
The following function implements the gradient magnitude using
Gaussian derivatives:

    The function :func:`gaussian_gradient_magnitude` calculates the
    gradient magnitude using :func:`gaussian_filter` to calculate the first
    derivatives. The standard-deviations of the Gaussian filter along
    each axis are passed through the parameter *sigma* as a sequence or
    numbers. If *sigma* is not a sequence but a single number, the
    standard deviation of the filter is equal along all directions.


.. _ndimage-genericfilters:

Generic filter functions
^^^^^^^^^^^^^^^^^^^^^^^^

To implement filter functions, generic functions can be used that accept a 
callable object that implements the filtering operation. The iteration over the 
input and output arrays is handled by these generic functions, along with such
details as the implementation of the boundary conditions. Only a
callable object implementing a callback function that does the
actual filtering work must be provided. The callback function can
also be written in C and passed using a :ctype:`PyCObject` (see
:ref:`ndimage-ccallbacks` for more information).

    The :func:`generic_filter1d` function implements a generic
    one-dimensional filter function, where the actual filtering
    operation must be supplied as a python function (or other callable
    object). The :func:`generic_filter1d` function iterates over the lines
    of an array and calls :func:`function` at each line. The arguments that
    are passed to :func:`function` are one-dimensional arrays of the
    :ctype:`tFloat64` type. The first contains the values of the current line.
    It is extended at the beginning end the end, according to the
    *filter_size* and *origin* arguments. The second array should be
    modified in-place to provide the output values of the line. For
    example consider a correlation along one dimension::

        >>> a = arange(12).reshape(3,4)
        >>> correlate1d(a, [1, 2, 3])
        array([[ 3,  8, 14, 17],
               [27, 32, 38, 41],
               [51, 56, 62, 65]])

    The same operation can be implemented using :func:`generic_filter1d` as
    follows::

        >>> def fnc(iline, oline):
        ...     oline[...] = iline[:-2] + 2 * iline[1:-1] + 3 * iline[2:]
        ... 
        >>> generic_filter1d(a, fnc, 3)
        array([[ 3,  8, 14, 17],
               [27, 32, 38, 41],
               [51, 56, 62, 65]])

    Here the origin of the kernel was (by default) assumed to be in the
    middle of the filter of length 3. Therefore, each input line was
    extended by one value at the beginning and at the end, before the
    function was called.

    Optionally extra arguments can be defined and passed to the filter
    function. The *extra_arguments* and *extra_keywords* arguments
    can be used to pass a tuple of extra arguments and/or a dictionary
    of named arguments that are passed to derivative at each call. For
    example, we can pass the parameters of our filter as an argument::

        >>> def fnc(iline, oline, a, b):
        ...     oline[...] = iline[:-2] + a * iline[1:-1] + b * iline[2:]
        ... 
        >>> generic_filter1d(a, fnc, 3, extra_arguments = (2, 3))
        array([[ 3,  8, 14, 17],
               [27, 32, 38, 41],
               [51, 56, 62, 65]])

    or::

        >>> generic_filter1d(a, fnc, 3, extra_keywords = {'a':2, 'b':3})
        array([[ 3,  8, 14, 17],
               [27, 32, 38, 41],
               [51, 56, 62, 65]])

    The :func:`generic_filter` function implements a generic filter
    function, where the actual filtering operation must be supplied as
    a python function (or other callable object). The :func:`generic_filter`
    function iterates over the array and calls :func:`function` at each
    element. The argument of :func:`function` is a one-dimensional array of
    the :ctype:`tFloat64` type, that contains the values around the current
    element that are within the footprint of the filter. The function
    should return a single value that can be converted to a double
    precision number. For example consider a correlation::

        >>> a = arange(12).reshape(3,4)
        >>> correlate(a, [[1, 0], [0, 3]])
        array([[ 0,  3,  7, 11],
               [12, 15, 19, 23],
               [28, 31, 35, 39]])

    The same operation can be implemented using *generic_filter* as
    follows::

        >>> def fnc(buffer): 
        ...     return (buffer * array([1, 3])).sum()
        ... 
        >>> generic_filter(a, fnc, footprint = [[1, 0], [0, 1]])
        array([[ 0  3  7 11],
               [12 15 19 23],
               [28 31 35 39]])

    Here a kernel footprint was specified that contains only two
    elements. Therefore the filter function receives a buffer of length
    equal to two, which was multiplied with the proper weights and the
    result summed.

    When calling :func:`generic_filter`, either the sizes of a rectangular
    kernel or the footprint of the kernel must be provided. The *size*
    parameter, if provided, must be a sequence of sizes or a single
    number in which case the size of the filter is assumed to be equal
    along each axis. The *footprint*, if provided, must be an array
    that defines the shape of the kernel by its non-zero elements.

    Optionally extra arguments can be defined and passed to the filter
    function. The *extra_arguments* and *extra_keywords* arguments
    can be used to pass a tuple of extra arguments and/or a dictionary
    of named arguments that are passed to derivative at each call. For
    example, we can pass the parameters of our filter as an argument::

        >>> def fnc(buffer, weights): 
        ...     weights = asarray(weights)
        ...     return (buffer * weights).sum()
        ... 
        >>> generic_filter(a, fnc, footprint = [[1, 0], [0, 1]], extra_arguments = ([1, 3],))
        array([[ 0,  3,  7, 11],
               [12, 15, 19, 23],
               [28, 31, 35, 39]])

    or::

        >>> generic_filter(a, fnc, footprint = [[1, 0], [0, 1]], extra_keywords= {'weights': [1, 3]})
        array([[ 0,  3,  7, 11],
               [12, 15, 19, 23],
               [28, 31, 35, 39]])


These functions iterate over the lines or elements starting at the
last axis, i.e. the last index changes the fastest. This order of
iteration is guaranteed for the case that it is important to adapt
the filter depending on spatial location. Here is an example of
using a class that implements the filter and keeps track of the
current coordinates while iterating. It performs the same filter
operation as described above for :func:`generic_filter`, but
additionally prints the current coordinates:

::

    >>> a = arange(12).reshape(3,4)
    >>> 
    >>> class fnc_class:
    ...     def __init__(self, shape):
    ...         # store the shape:
    ...         self.shape = shape
    ...         # initialize the coordinates:
    ...         self.coordinates = [0] * len(shape)
    ...         
    ...     def filter(self, buffer):
    ...         result = (buffer * array([1, 3])).sum()
    ...         print self.coordinates
    ...         # calculate the next coordinates:
    ...         axes = range(len(self.shape))
    ...         axes.reverse()
    ...         for jj in axes:
    ...             if self.coordinates[jj] < self.shape[jj] - 1:
    ...                 self.coordinates[jj] += 1
    ...                 break
    ...             else:
    ...                 self.coordinates[jj] = 0
    ...         return result
    ... 
    >>> fnc = fnc_class(shape = (3,4))
    >>> generic_filter(a, fnc.filter, footprint = [[1, 0], [0, 1]]) 
    [0, 0]
    [0, 1]
    [0, 2]
    [0, 3]
    [1, 0]
    [1, 1]
    [1, 2]
    [1, 3]
    [2, 0]
    [2, 1]
    [2, 2]
    [2, 3]
    array([[ 0,  3,  7, 11],
           [12, 15, 19, 23],
           [28, 31, 35, 39]])

For the :func:`generic_filter1d` function the same approach works,
except that this function does not iterate over the axis that is
being filtered. The example for :func:`generic_filter1d` then becomes
this:

::

    >>> a = arange(12).reshape(3,4)
    >>> 
    >>> class fnc1d_class:
    ...     def __init__(self, shape, axis = -1):
    ...         # store the filter axis:
    ...         self.axis = axis
    ...         # store the shape:
    ...         self.shape = shape
    ...         # initialize the coordinates:
    ...         self.coordinates = [0] * len(shape)
    ...         
    ...     def filter(self, iline, oline):
    ...         oline[...] = iline[:-2] + 2 * iline[1:-1] + 3 * iline[2:]
    ...         print self.coordinates
    ...         # calculate the next coordinates:
    ...         axes = range(len(self.shape))
    ...         # skip the filter axis:
    ...         del axes[self.axis]
    ...         axes.reverse()
    ...         for jj in axes:
    ...             if self.coordinates[jj] < self.shape[jj] - 1:
    ...                 self.coordinates[jj] += 1
    ...                 break
    ...             else:
    ...                 self.coordinates[jj] = 0
    ... 
    >>> fnc = fnc1d_class(shape = (3,4))
    >>> generic_filter1d(a, fnc.filter, 3)
    [0, 0]
    [1, 0]
    [2, 0]
    array([[ 3,  8, 14, 17],
           [27, 32, 38, 41],
           [51, 56, 62, 65]])

Fourier domain filters
^^^^^^^^^^^^^^^^^^^^^^

The functions described in this section perform filtering
operations in the Fourier domain. Thus, the input array of such a
function should be compatible with an inverse Fourier transform
function, such as the functions from the :mod:`numpy.fft` module. We
therefore have to deal with arrays that may be the result of a real
or a complex Fourier transform. In the case of a real Fourier
transform only half of the of the symmetric complex transform is
stored. Additionally, it needs to be known what the length of the
axis was that was transformed by the real fft. The functions
described here provide a parameter *n* that in the case of a real
transform must be equal to the length of the real transform axis
before transformation. If this parameter is less than zero, it is
assumed that the input array was the result of a complex Fourier
transform. The parameter *axis* can be used to indicate along which
axis the real transform was executed.

    The :func:`fourier_shift` function multiplies the input array with the
    multi-dimensional Fourier transform of a shift operation for the
    given shift. The *shift* parameter is a sequences of shifts for
    each dimension, or a single value for all dimensions.


    The :func:`fourier_gaussian` function multiplies the input array with
    the multi-dimensional Fourier transform of a Gaussian filter with
    given standard-deviations *sigma*. The *sigma* parameter is a
    sequences of values for each dimension, or a single value for all
    dimensions.


    The :func:`fourier_uniform` function multiplies the input array with the
    multi-dimensional Fourier transform of a uniform filter with given
    sizes *size*. The *size* parameter is a sequences of values for
    each dimension, or a single value for all dimensions.


    The :func:`fourier_ellipsoid` function multiplies the input array with
    the multi-dimensional Fourier transform of a elliptically shaped
    filter with given sizes *size*. The *size* parameter is a sequences
    of values for each dimension, or a single value for all dimensions.
    This function is only implemented for dimensions 1, 2, and 3.


.. _ndimage-interpolation:

Interpolation functions
-----------------------

.. currentmodule:: scipy.ndimage.interpolation


This section describes various interpolation functions that are
based on B-spline theory. A good introduction to B-splines can be
found in: M. Unser, "Splines: A Perfect Fit for Signal and Image
Processing," IEEE Signal Processing Magazine, vol. 16, no. 6, pp.
22-38, November 1999. 

Spline pre-filters
^^^^^^^^^^^^^^^^^^

Interpolation using
splines of an order larger than 1 requires a pre- filtering step.
The interpolation functions described in section
:ref:`ndimage-interpolation` apply pre-filtering by calling
:func:`spline_filter`, but they can be instructed not to do this by
setting the *prefilter* keyword equal to False. This is useful if
more than one interpolation operation is done on the same array. In
this case it is more efficient to do the pre-filtering only once
and use a prefiltered array as the input of the interpolation
functions. The following two functions implement the
pre-filtering:

    The :func:`spline_filter1d` function calculates a one-dimensional spline
    filter along the given axis. An output array can optionally be
    provided. The order of the spline must be larger then 1 and less
    than 6.


    The :func:`spline_filter` function calculates a multi-dimensional spline
    filter.

    .. note:: The multi-dimensional filter is implemented as a sequence of one-dimensional spline filters. The intermediate arrays are stored in the same data type as the output. Therefore, if an output with a limited precision is requested, the results may be imprecise because intermediate results may be stored with insufficient precision. This can be prevented by specifying a output type of high precision.


Interpolation functions
^^^^^^^^^^^^^^^^^^^^^^^

Following functions all employ spline interpolation to effect some type of 
geometric transformation of the input array. This requires a mapping of the 
output coordinates to the input coordinates, and therefore the possibility 
arises that input values outside the boundaries are needed. This problem is
solved in the same way as described in :ref:`ndimage-filter-functions` 
for the multi-dimensional filter functions. Therefore these functions all 
support a *mode* parameter that determines how the boundaries are handled, and 
a *cval* parameter that gives a constant value in case that the 'constant'
mode is used.

    The :func:`geometric_transform` function applies an arbitrary geometric
    transform to the input. The given *mapping* function is called at
    each point in the output to find the corresponding coordinates in
    the input. *mapping* must be a callable object that accepts a tuple
    of length equal to the output array rank and returns the
    corresponding input coordinates as a tuple of length equal to the
    input array rank. The output shape and output type can optionally
    be provided. If not given they are equal to the input shape and
    type.

    For example::

        >>> a = arange(12).reshape(4,3).astype(np.float64)
        >>> def shift_func(output_coordinates):
        ...     return (output_coordinates[0] - 0.5, output_coordinates[1] - 0.5)
        ... 
        >>> geometric_transform(a, shift_func)
        array([[ 0.    ,  0.    ,  0.    ],
               [ 0.    ,  1.3625,  2.7375],
               [ 0.    ,  4.8125,  6.1875],
               [ 0.    ,  8.2625,  9.6375]])

    Optionally extra arguments can be defined and passed to the filter
    function. The *extra_arguments* and *extra_keywords* arguments
    can be used to pass a tuple of extra arguments and/or a dictionary
    of named arguments that are passed to derivative at each call. For
    example, we can pass the shifts in our example as arguments::

        >>> def shift_func(output_coordinates, s0, s1):
        ...     return (output_coordinates[0] - s0, output_coordinates[1] - s1)
        ... 
        >>> geometric_transform(a, shift_func, extra_arguments = (0.5, 0.5))
        array([[ 0.    ,  0.    ,  0.    ],
               [ 0.    ,  1.3625,  2.7375],
               [ 0.    ,  4.8125,  6.1875],
               [ 0.    ,  8.2625,  9.6375]])

    or::

        >>> geometric_transform(a, shift_func, extra_keywords = {'s0': 0.5, 's1': 0.5})
        array([[ 0.    ,  0.    ,  0.    ],
               [ 0.    ,  1.3625,  2.7375],
               [ 0.    ,  4.8125,  6.1875],
               [ 0.    ,  8.2625,  9.6375]])

    .. note:: The mapping function can also be written in C and passed using a :ctype:`PyCObject`. See :ref:`ndimage-ccallbacks` for more information.


    The function :func:`map_coordinates` applies an arbitrary coordinate
    transformation using the given array of coordinates. The shape of
    the output is derived from that of the coordinate array by dropping
    the first axis. The parameter *coordinates* is used to find for
    each point in the output the corresponding coordinates in the
    input. The values of *coordinates* along the first axis are the
    coordinates in the input array at which the output value is found.
    (See also the numarray `coordinates` function.) Since the
    coordinates may be non- integer coordinates, the value of the input
    at these coordinates is determined by spline interpolation of the
    requested order. Here is an example that interpolates a 2D array at
    (0.5, 0.5) and (1, 2)::

        >>> a = arange(12).reshape(4,3).astype(np.float64)
        >>> a
        array([[  0.,   1.,   2.],
               [  3.,   4.,   5.],
               [  6.,   7.,   8.],
               [  9.,  10.,  11.]])
        >>> map_coordinates(a, [[0.5, 2], [0.5, 1]])
        array([ 1.3625  7.    ])

    The :func:`affine_transform` function applies an affine transformation
    to the input array. The given transformation *matrix* and *offset*
    are used to find for each point in the output the corresponding
    coordinates in the input. The value of the input at the calculated
    coordinates is determined by spline interpolation of the requested
    order. The transformation *matrix* must be two-dimensional or can
    also be given as a one-dimensional sequence or array. In the latter
    case, it is assumed that the matrix is diagonal. A more efficient
    interpolation algorithm is then applied that exploits the
    separability of the problem. The output shape and output type can
    optionally be provided. If not given they are equal to the input
    shape and type.

    The :func:`shift` function returns a shifted version of the input, using
    spline interpolation of the requested *order*.

    The :func:`zoom` function returns a rescaled version of the input, using
    spline interpolation of the requested *order*.

    The :func:`rotate` function returns the input array rotated in the plane
    defined by the two axes given by the parameter *axes*, using spline
    interpolation of the requested *order*. The angle must be given in
    degrees. If *reshape* is true, then the size of the output array is
    adapted to contain the rotated input.


.. _ndimage-morphology:

Morphology
----------

.. _ndimage-binary-morphology:

Binary morphology
^^^^^^^^^^^^^^^^^

.. currentmodule:: scipy.ndimage.morphology

Binary morphology (need something to put here).

    The :func:`generate_binary_structure` functions generates a binary
    structuring element for use in binary morphology operations. The
    *rank* of the structure must be provided. The size of the structure
    that is returned is equal to three in each direction. The value of
    each element is equal to one if the square of the Euclidean
    distance from the element to the center is less or equal to
    *connectivity*. For instance, two dimensional 4-connected and
    8-connected structures are generated as follows::

        >>> generate_binary_structure(2, 1)
        array([[False,  True, False],
               [ True,  True,  True],
               [False,  True, False]], dtype=bool)
        >>> generate_binary_structure(2, 2)
        array([[ True,  True,  True],
               [ True,  True,  True],
               [ True,  True,  True]], dtype=bool)

Most binary morphology functions can be expressed in terms of the
basic operations erosion and dilation:

    The :func:`binary_erosion` function implements binary erosion of arrays
    of arbitrary rank with the given structuring element. The origin
    parameter controls the placement of the structuring element as
    described in :ref:`ndimage-filter-functions`. If no
    structuring element is provided, an element with connectivity equal
    to one is generated using :func:`generate_binary_structure`. The
    *border_value* parameter gives the value of the array outside
    boundaries. The erosion is repeated *iterations* times. If
    *iterations* is less than one, the erosion is repeated until the
    result does not change anymore. If a *mask* array is given, only
    those elements with a true value at the corresponding mask element
    are modified at each iteration.

    The :func:`binary_dilation` function implements binary dilation of
    arrays of arbitrary rank with the given structuring element. The
    origin parameter controls the placement of the structuring element
    as described in :ref:`ndimage-filter-functions`. If no
    structuring element is provided, an element with connectivity equal
    to one is generated using :func:`generate_binary_structure`. The
    *border_value* parameter gives the value of the array outside
    boundaries. The dilation is repeated *iterations* times. If
    *iterations* is less than one, the dilation is repeated until the
    result does not change anymore. If a *mask* array is given, only
    those elements with a true value at the corresponding mask element
    are modified at each iteration.

    Here is an example of using :func:`binary_dilation` to find all elements
    that touch the border, by repeatedly dilating an empty array from
    the border using the data array as the mask::

        >>> struct = array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
        >>> a = array([[1,0,0,0,0], [1,1,0,1,0], [0,0,1,1,0], [0,0,0,0,0]])
        >>> a
        array([[1, 0, 0, 0, 0],
               [1, 1, 0, 1, 0],
               [0, 0, 1, 1, 0],
               [0, 0, 0, 0, 0]])
        >>> binary_dilation(zeros(a.shape), struct, -1, a, border_value=1)
        array([[ True, False, False, False, False],
               [ True,  True, False, False, False],
               [False, False, False, False, False],
               [False, False, False, False, False]], dtype=bool)


The :func:`binary_erosion` and :func:`binary_dilation` functions both have an
*iterations* parameter which allows the erosion or dilation to be
repeated a number of times. Repeating an erosion or a dilation with
a given structure *n* times is equivalent to an erosion or a
dilation with a structure that is *n-1* times dilated with itself.
A function is provided that allows the calculation of a structure
that is dilated a number of times with itself:

    The :func:`iterate_structure` function returns a structure by dilation
    of the input structure *iteration* - 1 times with itself. For
    instance::

        >>> struct = generate_binary_structure(2, 1)
        >>> struct
        array([[False,  True, False],
               [ True,  True,  True],
               [False,  True, False]], dtype=bool)
        >>> iterate_structure(struct, 2)
        array([[False, False,  True, False, False],
               [False,  True,  True,  True, False],
               [ True,  True,  True,  True,  True],
               [False,  True,  True,  True, False],
               [False, False,  True, False, False]], dtype=bool)

    If the origin of the original structure is equal to 0, then it is
    also equal to 0 for the iterated structure. If not, the origin must
    also be adapted if the equivalent of the *iterations* erosions or
    dilations must be achieved with the iterated structure. The adapted
    origin is simply obtained by multiplying with the number of
    iterations. For convenience the :func:`iterate_structure` also returns
    the adapted origin if the *origin* parameter is not None::

        >>> iterate_structure(struct, 2, -1)
        (array([[False, False,  True, False, False],
               [False,  True,  True,  True, False],
               [ True,  True,  True,  True,  True],
               [False,  True,  True,  True, False],
               [False, False,  True, False, False]], dtype=bool), [-2, -2])


Other morphology operations can be defined in terms of erosion and
d dilation. Following functions provide a few of these operations
for convenience:

    The :func:`binary_opening` function implements binary opening of arrays
    of arbitrary rank with the given structuring element. Binary
    opening is equivalent to a binary erosion followed by a binary
    dilation with the same structuring element. The origin parameter
    controls the placement of the structuring element as described in
    :ref:`ndimage-filter-functions`. If no structuring element is
    provided, an element with connectivity equal to one is generated
    using :func:`generate_binary_structure`. The *iterations* parameter
    gives the number of erosions that is performed followed by the same
    number of dilations.


    The :func:`binary_closing` function implements binary closing of arrays
    of arbitrary rank with the given structuring element. Binary
    closing is equivalent to a binary dilation followed by a binary
    erosion with the same structuring element. The origin parameter
    controls the placement of the structuring element as described in
    :ref:`ndimage-filter-functions`. If no structuring element is
    provided, an element with connectivity equal to one is generated
    using :func:`generate_binary_structure`. The *iterations* parameter
    gives the number of dilations that is performed followed by the
    same number of erosions.


    The :func:`binary_fill_holes` function is used to close holes in
    objects in a binary image, where the structure defines the
    connectivity of the holes. The origin parameter controls the
    placement of the structuring element as described in 
    :ref:`ndimage-filter-functions`. If no structuring element is
    provided, an element with connectivity equal to one is generated
    using :func:`generate_binary_structure`.


    The :func:`binary_hit_or_miss` function implements a binary
    hit-or-miss transform of arrays of arbitrary rank with the given
    structuring elements. The hit-or-miss transform is calculated by
    erosion of the input with the first structure, erosion of the
    logical *not* of the input with the second structure, followed by
    the logical *and* of these two erosions. The origin parameters
    control the placement of the structuring elements as described in
    :ref:`ndimage-filter-functions`. If *origin2* equals None it
    is set equal to the *origin1* parameter. If the first structuring
    element is not provided, a structuring element with connectivity
    equal to one is generated using :func:`generate_binary_structure`, if
    *structure2* is not provided, it is set equal to the logical *not*
    of *structure1*.


.. _ndimage-grey-morphology:

Grey-scale morphology
^^^^^^^^^^^^^^^^^^^^^

.. currentmodule:: scipy.ndimage.morphology

Grey-scale morphology operations are the equivalents of binary
morphology operations that operate on arrays with arbitrary values.
Below we describe the grey-scale equivalents of erosion, dilation,
opening and closing. These operations are implemented in a similar
fashion as the filters described in 
:ref:`ndimage-filter-functions`, and we refer to this section for the
description of filter kernels and footprints, and the handling of
array borders. The grey-scale morphology operations optionally take
a *structure* parameter that gives the values of the structuring
element. If this parameter is not given the structuring element is
assumed to be flat with a value equal to zero. The shape of the
structure can optionally be defined by the *footprint* parameter.
If this parameter is not given, the structure is assumed to be
rectangular, with sizes equal to the dimensions of the *structure*
array, or by the *size* parameter if *structure* is not given. The
*size* parameter is only used if both *structure* and *footprint*
are not given, in which case the structuring element is assumed to
be rectangular and flat with the dimensions given by *size*. The
*size* parameter, if provided, must be a sequence of sizes or a
single number in which case the size of the filter is assumed to be
equal along each axis. The *footprint* parameter, if provided, must
be an array that defines the shape of the kernel by its non-zero
elements.

Similar to binary erosion and dilation there are operations for
grey-scale erosion and dilation:

    The :func:`grey_erosion` function calculates a multi-dimensional grey-
    scale erosion.


    The :func:`grey_dilation` function calculates a multi-dimensional grey-
    scale dilation.


Grey-scale opening and closing operations can be defined similar to
their binary counterparts:

    The :func:`grey_opening` function implements grey-scale opening of
    arrays of arbitrary rank. Grey-scale opening is equivalent to a
    grey-scale erosion followed by a grey-scale dilation.


    The :func:`grey_closing` function implements grey-scale closing of
    arrays of arbitrary rank. Grey-scale opening is equivalent to a
    grey-scale dilation followed by a grey-scale erosion.


    The :func:`morphological_gradient` function implements a grey-scale
    morphological gradient of arrays of arbitrary rank. The grey-scale
    morphological gradient is equal to the difference of a grey-scale
    dilation and a grey-scale erosion.


    The :func:`morphological_laplace` function implements a grey-scale
    morphological laplace of arrays of arbitrary rank. The grey-scale
    morphological laplace is equal to the sum of a grey-scale dilation
    and a grey-scale erosion minus twice the input.


    The :func:`white_tophat` function implements a white top-hat filter of
    arrays of arbitrary rank. The white top-hat is equal to the
    difference of the input and a grey-scale opening.


    The :func:`black_tophat` function implements a black top-hat filter of
    arrays of arbitrary rank. The black top-hat is equal to the
    difference of the a grey-scale closing and the input.


.. _ndimage-distance-transforms:

Distance transforms
-------------------

.. currentmodule:: scipy.ndimage.morphology

Distance transforms are used to
calculate the minimum distance from each element of an object to
the background. The following functions implement distance
transforms for three different distance metrics: Euclidean, City
Block, and Chessboard distances.

    The function :func:`distance_transform_cdt` uses a chamfer type
    algorithm to calculate the distance transform of the input, by
    replacing each object element (defined by values larger than zero)
    with the shortest distance to the background (all non-object
    elements). The structure determines the type of chamfering that is
    done. If the structure is equal to 'cityblock' a structure is
    generated using :func:`generate_binary_structure` with a squared
    distance equal to 1. If the structure is equal to 'chessboard', a
    structure is generated using :func:`generate_binary_structure` with a
    squared distance equal to the rank of the array. These choices
    correspond to the common interpretations of the cityblock and the
    chessboard distancemetrics in two dimensions.

    In addition to the distance transform, the feature transform can be
    calculated. In this case the index of the closest background
    element is returned along the first axis of the result. The
    *return_distances*, and *return_indices* flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    The *distances* and *indices* arguments can be used to give
    optional output arrays that must be of the correct size and type
    (both :ctype:`Int32`).

    The basics of the algorithm used to implement this function is
    described in: G. Borgefors, "Distance transformations in arbitrary
    dimensions.", Computer Vision, Graphics, and Image Processing,
    27:321-345, 1984.


    The function :func:`distance_transform_edt` calculates the exact
    euclidean distance transform of the input, by replacing each object
    element (defined by values larger than zero) with the shortest
    euclidean distance to the background (all non-object elements).

    In addition to the distance transform, the feature transform can be
    calculated. In this case the index of the closest background
    element is returned along the first axis of the result. The
    *return_distances*, and *return_indices* flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    Optionally the sampling along each axis can be given by the
    *sampling* parameter which should be a sequence of length equal to
    the input rank, or a single number in which the sampling is assumed
    to be equal along all axes.

    The *distances* and *indices* arguments can be used to give
    optional output arrays that must be of the correct size and type
    (:ctype:`Float64` and :ctype:`Int32`).

    The algorithm used to implement this function is described in: C.
    R. Maurer, Jr., R. Qi, and V. Raghavan, "A linear time algorithm
    for computing exact euclidean distance transforms of binary images
    in arbitrary dimensions. IEEE Trans. PAMI 25, 265-270, 2003.


    The function :func:`distance_transform_bf` uses a brute-force algorithm
    to calculate the distance transform of the input, by replacing each
    object element (defined by values larger than zero) with the
    shortest distance to the background (all non-object elements). The
    metric must be one of "euclidean", "cityblock", or
    "chessboard".

    In addition to the distance transform, the feature transform can be
    calculated. In this case the index of the closest background
    element is returned along the first axis of the result. The
    *return_distances*, and *return_indices* flags can be used to
    indicate if the distance transform, the feature transform, or both
    must be returned.

    Optionally the sampling along each axis can be given by the
    *sampling* parameter which should be a sequence of length equal to
    the input rank, or a single number in which the sampling is assumed
    to be equal along all axes. This parameter is only used in the case
    of the euclidean distance transform.

    The *distances* and *indices* arguments can be used to give
    optional output arrays that must be of the correct size and type
    (:ctype:`Float64` and :ctype:`Int32`).

    .. note:: This function uses a slow brute-force algorithm, the function :func:`distance_transform_cdt` can be used to more efficiently calculate cityblock and chessboard distance transforms. The function :func:`distance_transform_edt` can be used to more efficiently calculate the exact euclidean distance transform.


Segmentation and labeling
-------------------------

Segmentation is the process of separating objects of interest from
the background. The most simple approach is probably intensity
thresholding, which is easily done with :mod:`numpy` functions::

    >>> a = array([[1,2,2,1,1,0],
    ...            [0,2,3,1,2,0],
    ...            [1,1,1,3,3,2],
    ...            [1,1,1,1,2,1]])
    >>> where(a > 1, 1, 0)
    array([[0, 1, 1, 0, 0, 0],
           [0, 1, 1, 0, 1, 0],
           [0, 0, 0, 1, 1, 1],
           [0, 0, 0, 0, 1, 0]])

The result is a binary image, in which the individual objects still
need to be identified and labeled. The function :func:`label` generates
an array where each object is assigned a unique number:

    The :func:`label` function generates an array where the objects in the
    input are labeled with an integer index. It returns a tuple
    consisting of the array of object labels and the number of objects
    found, unless the *output* parameter is given, in which case only
    the number of objects is returned. The connectivity of the objects
    is defined by a structuring element. For instance, in two
    dimensions using a four-connected structuring element gives::

        >>> a = array([[0,1,1,0,0,0],[0,1,1,0,1,0],[0,0,0,1,1,1],[0,0,0,0,1,0]])
        >>> s = [[0, 1, 0], [1,1,1], [0,1,0]]
        >>> label(a, s)
        (array([[0, 1, 1, 0, 0, 0],
               [0, 1, 1, 0, 2, 0],
               [0, 0, 0, 2, 2, 2],
               [0, 0, 0, 0, 2, 0]]), 2)

    These two objects are not connected because there is no way in
    which we can place the structuring element such that it overlaps
    with both objects. However, an 8-connected structuring element
    results in only a single object::

        >>> a = array([[0,1,1,0,0,0],[0,1,1,0,1,0],[0,0,0,1,1,1],[0,0,0,0,1,0]])
        >>> s = [[1,1,1], [1,1,1], [1,1,1]]
        >>> label(a, s)[0]
        array([[0, 1, 1, 0, 0, 0],
               [0, 1, 1, 0, 1, 0],
               [0, 0, 0, 1, 1, 1],
               [0, 0, 0, 0, 1, 0]])

    If no structuring element is provided, one is generated by calling
    :func:`generate_binary_structure` (see 
    :ref:`ndimage-binary-morphology`)
    using a connectivity of one (which in 2D is the 4-connected
    structure of the first example). The input can be of any type, any
    value not equal to zero is taken to be part of an object. This is
    useful if you need to 're-label' an array of object indices, for
    instance after removing unwanted objects. Just apply the label
    function again to the index array. For instance::

        >>> l, n = label([1, 0, 1, 0, 1])
        >>> l
        array([1 0 2 0 3])
        >>> l = where(l != 2, l, 0)
        >>> l
        array([1 0 0 0 3])
        >>> label(l)[0]
        array([1 0 0 0 2])

    .. note:: The structuring element used by :func:`label` is assumed to be symmetric.


There is a large number of other approaches for segmentation, for
instance from an estimation of the borders of the objects that can
be obtained for instance by derivative filters. One such an
approach is watershed segmentation. The function :func:`watershed_ift`
generates an array where each object is assigned a unique label,
from an array that localizes the object borders, generated for
instance by a gradient magnitude filter. It uses an array
containing initial markers for the objects:

    The :func:`watershed_ift` function applies a watershed from markers
    algorithm, using an Iterative Forest Transform, as described in: P.
    Felkel, R. Wegenkittl, and M. Bruckschwaiger, "Implementation and
    Complexity of the Watershed-from-Markers Algorithm Computed as a
    Minimal Cost Forest.", Eurographics 2001, pp. C:26-35.

    The inputs of this function are the array to which the transform is
    applied, and an array of markers that designate the objects by a
    unique label, where any non-zero value is a marker. For instance::

        >>> input = array([[0, 0, 0, 0, 0, 0, 0],
        ...                [0, 1, 1, 1, 1, 1, 0],
        ...                [0, 1, 0, 0, 0, 1, 0],
        ...                [0, 1, 0, 0, 0, 1, 0],
        ...                [0, 1, 0, 0, 0, 1, 0],
        ...                [0, 1, 1, 1, 1, 1, 0],
        ...                [0, 0, 0, 0, 0, 0, 0]], np.uint8)
        >>> markers = array([[1, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 2, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0]], np.int8)
        >>> watershed_ift(input, markers)
        array([[1, 1, 1, 1, 1, 1, 1],
               [1, 1, 2, 2, 2, 1, 1],
               [1, 2, 2, 2, 2, 2, 1],
               [1, 2, 2, 2, 2, 2, 1],
               [1, 2, 2, 2, 2, 2, 1],
               [1, 1, 2, 2, 2, 1, 1],
               [1, 1, 1, 1, 1, 1, 1]], dtype=int8)

    Here two markers were used to designate an object (*marker* = 2) and
    the background (*marker* = 1). The order in which these are processed
    is arbitrary: moving the marker for the background to the lower
    right corner of the array yields a different result::

        >>> markers = array([[0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 2, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 1]], np.int8)
        >>> watershed_ift(input, markers)
        array([[1, 1, 1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1, 1, 1],
               [1, 1, 2, 2, 2, 1, 1],
               [1, 1, 2, 2, 2, 1, 1],
               [1, 1, 2, 2, 2, 1, 1],
               [1, 1, 1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1, 1, 1]], dtype=int8)

    The result is that the object (*marker* = 2) is smaller because the
    second marker was processed earlier. This may not be the desired
    effect if the first marker was supposed to designate a background
    object. Therefore :func:`watershed_ift` treats markers with a negative
    value explicitly as background markers and processes them after the
    normal markers. For instance, replacing the first marker by a
    negative marker gives a result similar to the first example::

        >>> markers = array([[0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 2, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, 0],
        ...                  [0, 0, 0, 0, 0, 0, -1]], np.int8)
        >>> watershed_ift(input, markers)
        array([[-1, -1, -1, -1, -1, -1, -1],
               [-1, -1,  2,  2,  2, -1, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1, -1,  2,  2,  2, -1, -1],
               [-1, -1, -1, -1, -1, -1, -1]], dtype=int8)

    The connectivity of the objects is defined by a structuring
    element. If no structuring element is provided, one is generated by
    calling :func:`generate_binary_structure` (see
    :ref:`ndimage-binary-morphology`) using a connectivity of one 
    (which in 2D is a 4-connected structure.) For example, using 
    an 8-connected structure with the last example yields a different object::

        >>> watershed_ift(input, markers,
        ...               structure = [[1,1,1], [1,1,1], [1,1,1]])
        array([[-1, -1, -1, -1, -1, -1, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1,  2,  2,  2,  2,  2, -1],
               [-1, -1, -1, -1, -1, -1, -1]], dtype=int8)

    .. note:: The implementation of :func:`watershed_ift` limits the data types of the input to :ctype:`UInt8` and :ctype:`UInt16`.


.. _ndimage-object-measurements:

Object measurements
-------------------

.. currentmodule:: scipy.ndimage.measurements

Given an array of labeled objects, the properties of the individual
objects can be measured. The :func:`find_objects` function can be used
to generate a list of slices that for each object, give the
smallest sub-array that fully contains the object:

    The :func:`find_objects` function finds all objects in a labeled array and
    returns a list of slices that correspond to the smallest regions in
    the array that contains the object. For instance::

        >>> a = array([[0,1,1,0,0,0],[0,1,1,0,1,0],[0,0,0,1,1,1],[0,0,0,0,1,0]])
        >>> l, n = label(a)
        >>> f = find_objects(l)
        >>> a[f[0]]
        array([[1 1],
               [1 1]])
        >>> a[f[1]]
        array([[0, 1, 0],
               [1, 1, 1],
               [0, 1, 0]])

    :func:`find_objects` returns slices for all objects, unless the
    *max_label* parameter is larger then zero, in which case only the
    first *max_label* objects are returned. If an index is missing in
    the *label* array, None is return instead of a slice. For
    example::

        >>> find_objects([1, 0, 3, 4], max_label = 3)
        [(slice(0, 1, None),), None, (slice(2, 3, None),)]


The list of slices generated by :func:`find_objects` is useful to find
the position and dimensions of the objects in the array, but can
also be used to perform measurements on the individual objects. Say
we want to find the sum of the intensities of an object in image::

    >>> image = arange(4 * 6).reshape(4, 6)
    >>> mask = array([[0,1,1,0,0,0],[0,1,1,0,1,0],[0,0,0,1,1,1],[0,0,0,0,1,0]])
    >>> labels = label(mask)[0]
    >>> slices = find_objects(labels)

Then we can calculate the sum of the elements in the second
object::

    >>> where(labels[slices[1]] == 2, image[slices[1]], 0).sum()
    80

That is however not particularly efficient, and may also be more
complicated for other types of measurements. Therefore a few
measurements functions are defined that accept the array of object
labels and the index of the object to be measured. For instance
calculating the sum of the intensities can be done by::

    >>> sum(image, labels, 2)
    80

For large arrays and small objects it is more efficient to call the
measurement functions after slicing the array::

    >>> sum(image[slices[1]], labels[slices[1]], 2)
    80

Alternatively, we can do the measurements for a number of labels
with a single function call, returning a list of results. For
instance, to measure the sum of the values of the background and
the second object in our example we give a list of labels::

    >>> sum(image, labels, [0, 2])
    array([178.0, 80.0])

The measurement functions described below all support the *index*
parameter to indicate which object(s) should be measured. The
default value of *index* is None. This indicates that all
elements where the label is larger than zero should be treated as a
single object and measured. Thus, in this case the *labels* array
is treated as a mask defined by the elements that are larger than
zero. If *index* is a number or a sequence of numbers it gives the
labels of the objects that are measured. If *index* is a sequence,
a list of the results is returned. Functions that return more than
one result, return their result as a tuple if *index* is a single
number, or as a tuple of lists, if *index* is a sequence.

    The :func:`sum` function calculates the sum of the elements of the object
    with label(s) given by *index*, using the *labels* array for the
    object labels. If *index* is None, all elements with a non-zero
    label value are treated as a single object. If *label* is None,
    all elements of *input* are used in the calculation.


    The :func:`mean` function calculates the mean of the elements of the
    object with label(s) given by *index*, using the *labels* array for
    the object labels. If *index* is None, all elements with a
    non-zero label value are treated as a single object. If *label* is
    None, all elements of *input* are used in the calculation.


    The :func:`variance` function calculates the variance of the elements of
    the object with label(s) given by *index*, using the *labels* array
    for the object labels. If *index* is None, all elements with a
    non-zero label value are treated as a single object. If *label* is
    None, all elements of *input* are used in the calculation.


    The :func:`standard_deviation` function calculates the standard
    deviation of the elements of the object with label(s) given by
    *index*, using the *labels* array for the object labels. If *index*
    is None, all elements with a non-zero label value are treated as
    a single object. If *label* is None, all elements of *input* are
    used in the calculation.


    The :func:`minimum` function calculates the minimum of the elements of
    the object with label(s) given by *index*, using the *labels* array
    for the object labels. If *index* is None, all elements with a
    non-zero label value are treated as a single object. If *label* is
    None, all elements of *input* are used in the calculation.


    The :func:`maximum` function calculates the maximum of the elements of
    the object with label(s) given by *index*, using the *labels* array
    for the object labels. If *index* is None, all elements with a
    non-zero label value are treated as a single object. If *label* is
    None, all elements of *input* are used in the calculation.


    The :func:`minimum_position` function calculates the position of the
    minimum of the elements of the object with label(s) given by
    *index*, using the *labels* array for the object labels. If *index*
    is None, all elements with a non-zero label value are treated as
    a single object. If *label* is None, all elements of *input* are
    used in the calculation.


    The :func:`maximum_position` function calculates the position of the
    maximum of the elements of the object with label(s) given by
    *index*, using the *labels* array for the object labels. If *index*
    is None, all elements with a non-zero label value are treated as
    a single object. If *label* is None, all elements of *input* are
    used in the calculation.


    The :func:`extrema` function calculates the minimum, the maximum, and
    their positions, of the elements of the object with label(s) given
    by *index*, using the *labels* array for the object labels. If
    *index* is None, all elements with a non-zero label value are
    treated as a single object. If *label* is None, all elements of
    *input* are used in the calculation. The result is a tuple giving
    the minimum, the maximum, the position of the minimum and the
    postition of the maximum. The result is the same as a tuple formed
    by the results of the functions *minimum*, *maximum*,
    *minimum_position*, and *maximum_position* that are described
    above.


    The :func:`center_of_mass` function calculates the center of mass of
    the of the object with label(s) given by *index*, using the
    *labels* array for the object labels. If *index* is None, all
    elements with a non-zero label value are treated as a single
    object. If *label* is None, all elements of *input* are used in
    the calculation.


    The :func:`histogram` function calculates a histogram of the of the
    object with label(s) given by *index*, using the *labels* array for
    the object labels. If *index* is None, all elements with a
    non-zero label value are treated as a single object. If *label* is
    None, all elements of *input* are used in the calculation.
    Histograms are defined by their minimum (*min*), maximum (*max*)
    and the number of bins (*bins*). They are returned as
    one-dimensional arrays of type :ctype:`Int32`.


.. _ndimage-ccallbacks:

Extending :mod:`ndimage` in C
-----------------------------

.. highlight:: c

A few functions in the :mod:`scipy.ndimage` take a call-back 
argument. This can be a python function, but also a :ctype:`PyCObject`
containing a pointer to a C function. To use this feature, you must 
write your own C extension that defines the function, and define a Python function that returns a :ctype:`PyCObject` containing a pointer to this function.

An example of a function that supports this is
:func:`geometric_transform` (see :ref:`ndimage-interpolation`).
You can pass it a python callable object that defines a mapping
from all output coordinates to corresponding coordinates in the
input array. This mapping function can also be a C function, which
generally will be much more efficient, since the overhead of
calling a python function at each element is avoided.

For example to implement a simple shift function we define the
following function:

::

    static int 
    _shift_function(int *output_coordinates, double* input_coordinates,
                    int output_rank, int input_rank, void *callback_data)
    {
      int ii;
      /* get the shift from the callback data pointer: */
      double shift = *(double*)callback_data;
      /* calculate the coordinates: */
      for(ii = 0; ii < irank; ii++)
        icoor[ii] = ocoor[ii] - shift;
      /* return OK status: */
      return 1;
    }

This function is called at every element of the output array,
passing the current coordinates in the *output_coordinates* array.
On return, the *input_coordinates* array must contain the
coordinates at which the input is interpolated. The ranks of the
input and output array are passed through *output_rank* and
*input_rank*. The value of the shift is passed through the
*callback_data* argument, which is a pointer to void. The function
returns an error status, in this case always 1, since no error can
occur.

A pointer to this function and a pointer to the shift value must be
passed to :func:`geometric_transform`. Both are passed by a single
:ctype:`PyCObject` which is created by the following python extension
function:

::

    static PyObject *
    py_shift_function(PyObject *obj, PyObject *args)
    {
      double shift = 0.0;
      if (!PyArg_ParseTuple(args, "d", &shift)) {
        PyErr_SetString(PyExc_RuntimeError, "invalid parameters");
        return NULL;
      } else {
        /* assign the shift to a dynamically allocated location: */
        double *cdata = (double*)malloc(sizeof(double));
        *cdata = shift;
        /* wrap function and callback_data in a CObject: */
        return PyCObject_FromVoidPtrAndDesc(_shift_function, cdata,
                                            _destructor);
      }
    }

The value of the shift is obtained and then assigned to a
dynamically allocated memory location. Both this data pointer and
the function pointer are then wrapped in a :ctype:`PyCObject`, which is
returned. Additionally, a pointer to a destructor function is
given, that will free the memory we allocated for the shift value
when the :ctype:`PyCObject` is destroyed. This destructor is very simple:

::

    static void
    _destructor(void* cobject, void *cdata)
    {
      if (cdata)
        free(cdata);
    }

To use these functions, an extension module is built:

::

    static PyMethodDef methods[] = {
      {"shift_function", (PyCFunction)py_shift_function, METH_VARARGS, ""},
      {NULL, NULL, 0, NULL}
    };
    
    void
    initexample(void)
    {
      Py_InitModule("example", methods);
    }

This extension can then be used in Python, for example:

.. highlight:: python

::

    >>> import example
    >>> array = arange(12).reshape=(4, 3).astype(np.float64)
    >>> fnc = example.shift_function(0.5)
    >>> geometric_transform(array, fnc)
    array([[ 0.      0.      0.    ],
           [ 0.      1.3625  2.7375],
           [ 0.      4.8125  6.1875],
           [ 0.      8.2625  9.6375]])

C callback functions for use with :mod:`ndimage` functions must all
be written according to this scheme. The next section lists the
:mod:`ndimage` functions that acccept a C callback function and
gives the prototype of the callback function.

Functions that support C callback functions
-------------------------------------------

The :mod:`ndimage` functions that support C callback functions are
described here. Obviously, the prototype of the function that is
provided to these functions must match exactly that what they
expect. Therefore we give here the prototypes of the callback
functions. All these callback functions accept a void
*callback_data* pointer that must be wrapped in a :ctype:`PyCObject` using
the Python :cfunc:`PyCObject_FromVoidPtrAndDesc` function, which can also
accept a pointer to a destructor function to free any memory
allocated for *callback_data*. If *callback_data* is not needed,
:cfunc:`PyCObject_FromVoidPtr` may be used instead. The callback
functions must return an integer error status that is equal to zero
if something went wrong, or 1 otherwise. If an error occurs, you
should normally set the python error status with an informative
message before returning, otherwise, a default error message is set
by the calling function.

The function :func:`generic_filter` (see
:ref:`ndimage-genericfilters`) accepts a callback function with the
following prototype:

    The calling function iterates over the elements of the input and
    output arrays, calling the callback function at each element. The
    elements within the footprint of the filter at the current element
    are passed through the *buffer* parameter, and the number of
    elements within the footprint through *filter_size*. The
    calculated valued should be returned in the *return_value*
    argument.

The function :func:`generic_filter1d` (see
:ref:`ndimage-genericfilters`) accepts a callback function with the
following prototype:

    The calling function iterates over the lines of the input and
    output arrays, calling the callback function at each line. The
    current line is extended according to the border conditions set by
    the calling function, and the result is copied into the array that
    is passed through the *input_line* array. The length of the input
    line (after extension) is passed through *input_length*. The
    callback function should apply the 1D filter and store the result
    in the array passed through *output_line*. The length of the
    output line is passed through *output_length*.

The function :func:`geometric_transform` (see
:ref:`ndimage-interpolation`) expects a function with the following
prototype:

    The calling function iterates over the elements of the output
    array, calling the callback function at each element. The
    coordinates of the current output element are passed through
    *output_coordinates*. The callback function must return the
    coordinates at which the input must be interpolated in
    *input_coordinates*. The rank of the input and output arrays are
    given by *input_rank* and *output_rank* respectively.



