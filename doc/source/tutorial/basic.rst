Basic functions
===============

.. sectionauthor:: Travis E. Oliphant

.. currentmodule:: numpy

.. contents::

Interaction with NumPy
----------------------

SciPy builds on NumPy, and for all basic array handling needs you can
use NumPy functions:

    >>> import numpy as np
    >>> np.some_function()

Rather than giving a detailed description of each of these functions
(which is available in the NumPy Reference Guide or by using the
:func:`help`, :func:`info` and :func:`source` commands), this tutorial
will discuss some of the more useful commands, which require a little
introduction to use to their full potential.

To use functions from some of the SciPy modules, you can do:

    >>> from scipy import some_module
    >>> some_module.some_function()

The top level of :mod:`scipy` also contains functions from
:mod:`numpy` and :mod:`numpy.lib.scimath`. However, it is better to
use them directly from the :mod:`numpy` module instead.

Index tricks
^^^^^^^^^^^^

There are some class instances that make special use of the slicing
functionality to provide efficient means for array construction. This
part will discuss the operation of :obj:`numpy.mgrid` , :obj:`numpy.ogrid` ,
:obj:`numpy.r_` , and :obj:`numpy.c_` for quickly constructing arrays.

For example, rather than writing something like the following

    >>> a = np.concatenate(([3], [0]*5, np.arange(-1, 1.002, 2/9.0)))

with the :obj:`r_` command one can enter this as

    >>> a = np.r_[3,[0]*5,-1:1:10j]

which can ease typing and make for more readable code. Notice how
objects are concatenated, and the slicing syntax is (ab)used to
construct ranges. The other term that deserves a little explanation is
the use of the complex number 10j as the step size in the slicing
syntax. This non-standard use allows the number to be interpreted as
the number of points to produce in the range rather than as a step
size (note we would have used the long integer notation, 10L, but this
notation may go away in Python as the integers become unified). This
non-standard usage may be unsightly to some, but it gives the user the
ability to quickly construct complicated vectors in a very readable
fashion. When the number of points is specified in this way, the end-
point is inclusive.

The "r" stands for row concatenation because if the objects between
commas are 2-D arrays, they are stacked by rows (and thus
must have commensurate columns). There is an equivalent command
:obj:`c_` that stacks 2-D arrays by columns but works identically to
:obj:`r_` for 1-D arrays.

Another very useful class instance which makes use of extended slicing
notation is the function :obj:`mgrid`. In the simplest case, this
function can be used to construct 1-D ranges as a convenient substitute
for arange. It also allows the use of complex-numbers in the step-size
to indicate the number of points to place between the (inclusive)
end-points. The real purpose of this function however is to produce N,
N-D arrays, which provide coordinate arrays for an N-D
volume. The easiest way to understand this is with an example of its
usage:

    >>> np.mgrid[0:5,0:5]
    array([[[0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2],
            [3, 3, 3, 3, 3],
            [4, 4, 4, 4, 4]],
           [[0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4]]])
    >>> np.mgrid[0:5:4j,0:5:4j]
    array([[[ 0.    ,  0.    ,  0.    ,  0.    ],
            [ 1.6667,  1.6667,  1.6667,  1.6667],
            [ 3.3333,  3.3333,  3.3333,  3.3333],
            [ 5.    ,  5.    ,  5.    ,  5.    ]],
           [[ 0.    ,  1.6667,  3.3333,  5.    ],
            [ 0.    ,  1.6667,  3.3333,  5.    ],
            [ 0.    ,  1.6667,  3.3333,  5.    ],
            [ 0.    ,  1.6667,  3.3333,  5.    ]]])

Having meshed arrays like this is sometimes very useful. However, it
is not always needed just to evaluate some N-D function over
a grid due to the array-broadcasting rules of NumPy and SciPy. If this
is the only purpose for generating a meshgrid, you should instead use
the function :obj:`ogrid` which generates an "open" grid using :obj:`newaxis`
judiciously to create N, N-D arrays where only one dimension in each
array has length greater than 1. This will save memory and create the
same result if the only purpose for the meshgrid is to generate sample
points for evaluation of an N-D function.


Shape manipulation
^^^^^^^^^^^^^^^^^^

In this category of functions are routines for squeezing out length-
one dimensions from N-D arrays, ensuring that an array is at
least 1-, 2-, or 3-D, and stacking (concatenating) arrays by
rows, columns, and "pages" (in the third dimension). Routines for
splitting arrays (roughly the opposite of stacking arrays) are also
available.


Polynomials
^^^^^^^^^^^

There are two (interchangeable) ways to deal with 1-D polynomials in
SciPy. The first is to use the :class:`poly1d` class from NumPy. This
class accepts coefficients or polynomial roots to initialize a
polynomial. The polynomial object can then be manipulated in algebraic
expressions, integrated, differentiated, and evaluated. It even prints
like a polynomial:

    >>> from numpy import poly1d
    >>> p = poly1d([3,4,5])
    >>> print(p)
       2
    3 x + 4 x + 5
    >>> print(p*p)
       4      3      2
    9 x + 24 x + 46 x + 40 x + 25
    >>> print(p.integ(k=6))
       3     2
    1 x + 2 x + 5 x + 6
    >>> print(p.deriv())
    6 x + 4
    >>> p([4, 5])
    array([ 69, 100])

The other way to handle polynomials is as an array of coefficients
with the first element of the array giving the coefficient of the
highest power. There are explicit functions to add, subtract,
multiply, divide, integrate, differentiate, and evaluate polynomials
represented as sequences of coefficients.


Vectorizing functions (vectorize)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the features that NumPy provides is a class :obj:`vectorize` to
convert an ordinary Python function which accepts scalars and returns
scalars into a "vectorized-function" with the same broadcasting rules
as other NumPy functions (i.e., the Universal functions, or
ufuncs). For example, suppose you have a Python function named
``addsubtract`` defined as:

    >>> def addsubtract(a,b):
    ...    if a > b:
    ...        return a - b
    ...    else:
    ...        return a + b

which defines a function of two scalar variables and returns a scalar
result. The class vectorize can be used to "vectorize" this function so that ::

    >>> vec_addsubtract = np.vectorize(addsubtract)

returns a function which takes array arguments and returns an array
result:

    >>> vec_addsubtract([0,3,6,9],[1,3,5,7])
    array([1, 6, 1, 2])

This particular function could have been written in vector form
without the use of :obj:`vectorize`. However, functions that employ optimization
or integration routines can likely only be vectorized using ``vectorize.``

Type handling
^^^^^^^^^^^^^

Note the difference between :func:`numpy.iscomplex`/:func:`numpy.isreal` and
:func:`numpy.iscomplexobj`/:func:`numpy.isrealobj`. The former command is
array-based and returns byte arrays of ones and zeros providing the
result of the element-wise test. The latter command is object-based
and returns a scalar describing the result of the test on the entire
object.

Often it is required to get just the real and/or imaginary part of a
complex number. While complex numbers and arrays have attributes that
return those values, if one is not sure whether or not the object will
be complex-valued, it is better to use the functional forms
:func:`numpy.real` and :func:`numpy.imag` . These functions succeed for anything
that can be turned into a NumPy array. Consider also the function
:func:`numpy.real_if_close` which transforms a complex-valued number with a
tiny imaginary part into a real number.

Occasionally the need to check whether or not a number is a scalar
(Python (long)int, Python float, Python complex, or rank-0 array)
occurs in coding. This functionality is provided in the convenient
function :func:`numpy.isscalar` which returns a 1 or a 0.

Other useful functions
^^^^^^^^^^^^^^^^^^^^^^

There are also several other useful functions which should be
mentioned. For doing phase processing, the functions :func:`angle`,
and :obj:`unwrap` are useful. Also, the :obj:`linspace` and
:obj:`logspace` functions return equally spaced samples in a linear or
log scale.  Finally, it's useful to be aware of the indexing
capabilities of NumPy. Mention should be made of the function
:obj:`~numpy.select` which extends the functionality of :obj:`where` to
include multiple conditions and multiple choices. The calling
convention is ``select(condlist, choicelist, default=0)``. :obj:`numpy.select`
is a vectorized form of the multiple if-statement. It allows rapid
construction of a function which returns an array of results based on
a list of conditions. Each element of the return array is taken from
the array in a ``choicelist`` corresponding to the first condition in
``condlist`` that is true. For example:

    >>> x = np.arange(10)
    >>> condlist = [x<3, x>5]
    >>> choicelist = [x, x**2]
    >>> np.select(condlist, choicelist)
    array([ 0,  1,  2,  0,  0,  0, 36, 49, 64, 81])

Some additional useful functions can also be found in the module
:mod:`scipy.special`. For example the :obj:`~scipy.special.factorial`
and :obj:`~scipy.special.comb`
functions compute :math:`n!` and :math:`n!/k!(n-k)!` using either
exact integer arithmetic (thanks to Python's Long integer object), or
by using floating-point precision and the gamma function.

.. currentmodule:: scipy.misc

Other useful functions can be found in :mod:`scipy.misc`.
For example, two functions are provided that are useful for approximating
derivatives of functions using discrete-differences. The function
:obj:`central_diff_weights` returns weighting coefficients for an
equally-spaced :math:`N`-point approximation to the derivative of
order *o*. These weights must be multiplied by the function
corresponding to these points and the results added to obtain the
derivative approximation. This function is intended for use when only
samples of the function are available. When the function is an object
that can be handed to a routine and evaluated, the function
:obj:`derivative` can be used to automatically evaluate the object at
the correct points to obtain an N-point approximation to the *o*-th
derivative at a given point.
