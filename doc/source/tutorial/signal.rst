Signal Processing (`scipy.signal`)
==================================

.. sectionauthor:: Travis E. Oliphant

.. sectionauthor:: Pim Schellart

.. currentmodule:: scipy.signal

The signal processing toolbox currently contains some filtering
functions, a limited set of filter design tools, and a few B-spline
interpolation algorithms for one- and two-dimensional data. While the
B-spline algorithms could technically be placed under the
interpolation category, they are included here because they only work
with equally-spaced data and make heavy use of filter-theory and
transfer-function formalism to provide a fast B-spline transform. To
understand this section you will need to understand that a signal in
SciPy is an array of real or complex numbers.


B-splines
---------

A B-spline is an approximation of a continuous function over a finite-
domain in terms of B-spline coefficients and knot points. If the knot-
points are equally spaced with spacing :math:`\Delta x` , then the B-spline
approximation to a 1-dimensional function is the finite-basis expansion.

.. math::
   :nowrap:

    \[ y\left(x\right)\approx\sum_{j}c_{j}\beta^{o}\left(\frac{x}{\Delta x}-j\right).\]

In two dimensions with knot-spacing :math:`\Delta x` and :math:`\Delta y` , the
function representation is

.. math::
   :nowrap:

    \[ z\left(x,y\right)\approx\sum_{j}\sum_{k}c_{jk}\beta^{o}\left(\frac{x}{\Delta x}-j\right)\beta^{o}\left(\frac{y}{\Delta y}-k\right).\]

In these expressions, :math:`\beta^{o}\left(\cdot\right)` is the space-limited
B-spline basis function of order, :math:`o` . The requirement of equally-spaced
knot-points and equally-spaced data points, allows the development of fast
(inverse-filtering) algorithms for determining the coefficients, :math:`c_{j}`
, from sample-values, :math:`y_{n}` . Unlike the general spline interpolation
algorithms, these algorithms can quickly find the spline coefficients for large
images.

The advantage of representing a set of samples via B-spline basis
functions is that continuous-domain operators (derivatives, re-
sampling, integral, etc.) which assume that the data samples are drawn
from an underlying continuous function can be computed with relative
ease from the spline coefficients. For example, the second-derivative
of a spline is

.. math::
   :nowrap:

    \[ y{}^{\prime\prime}\left(x\right)=\frac{1}{\Delta x^{2}}\sum_{j}c_{j}\beta^{o\prime\prime}\left(\frac{x}{\Delta x}-j\right).\]

Using the property of B-splines that

.. math::
   :nowrap:

    \[ \frac{d^{2}\beta^{o}\left(w\right)}{dw^{2}}=\beta^{o-2}\left(w+1\right)-2\beta^{o-2}\left(w\right)+\beta^{o-2}\left(w-1\right)\]

it can be seen that

.. math::
   :nowrap:

    \[ y^{\prime\prime}\left(x\right)=\frac{1}{\Delta x^{2}}\sum_{j}c_{j}\left[\beta^{o-2}\left(\frac{x}{\Delta x}-j+1\right)-2\beta^{o-2}\left(\frac{x}{\Delta x}-j\right)+\beta^{o-2}\left(\frac{x}{\Delta x}-j-1\right)\right].\]

If :math:`o=3` , then at the sample points,

.. math::
   :nowrap:

    \begin{eqnarray*} \Delta x^{2}\left.y^{\prime}\left(x\right)\right|_{x=n\Delta x} & = & \sum_{j}c_{j}\delta_{n-j+1}-2c_{j}\delta_{n-j}+c_{j}\delta_{n-j-1},\\  & = & c_{n+1}-2c_{n}+c_{n-1}.\end{eqnarray*}

Thus, the second-derivative signal can be easily calculated from the
spline fit. if desired, smoothing splines can be found to make the
second-derivative less sensitive to random-errors.

The savvy reader will have already noticed that the data samples are
related to the knot coefficients via a convolution operator, so that
simple convolution with the sampled B-spline function recovers the
original data from the spline coefficients. The output of convolutions
can change depending on how boundaries are handled (this becomes
increasingly more important as the number of dimensions in the data-
set increases). The algorithms relating to B-splines in the signal-
processing sub package assume mirror-symmetric boundary conditions.
Thus, spline coefficients are computed based on that assumption, and
data-samples can be recovered exactly from the spline coefficients by
assuming them to be mirror-symmetric also.

Currently the package provides functions for determining second- and
third-order cubic spline coefficients from equally spaced samples in
one- and two-dimensions (:func:`signal.qspline1d`,
:func:`signal.qspline2d`, :func:`signal.cspline1d`,
:func:`signal.cspline2d`). The package also supplies a function (
:obj:`signal.bspline` ) for evaluating the bspline basis function,
:math:`\beta^{o}\left(x\right)` for arbitrary order and :math:`x.` For
large :math:`o` , the B-spline basis function can be approximated well
by a zero-mean Gaussian function with standard-deviation equal to
:math:`\sigma_{o}=\left(o+1\right)/12` :

.. math::
   :nowrap:

    \[ \beta^{o}\left(x\right)\approx\frac{1}{\sqrt{2\pi\sigma_{o}^{2}}}\exp\left(-\frac{x^{2}}{2\sigma_{o}}\right).\]

A function to compute this Gaussian for arbitrary :math:`x` and
:math:`o` is also available ( :obj:`signal.gauss_spline` ). The
following code and Figure uses spline-filtering to compute an
edge-image (the second-derivative of a smoothed spline) of Lena's face
which is an array returned by the command :func:`lena`. The command
:obj:`signal.sepfir2d` was used to apply a separable two-dimensional
FIR filter with mirror- symmetric boundary conditions to the spline
coefficients. This function is ideally suited for reconstructing
samples from spline coefficients and is faster than
:obj:`signal.convolve2d` which convolves arbitrary two-dimensional
filters and allows for choosing mirror-symmetric boundary conditions.

.. plot::

   >>> from numpy import *
   >>> from scipy import signal, misc
   >>> import matplotlib.pyplot as plt

   >>> image = misc.lena().astype(float32)
   >>> derfilt = array([1.0,-2,1.0],float32)
   >>> ck = signal.cspline2d(image,8.0)
   >>> deriv = signal.sepfir2d(ck, derfilt, [1]) + \
   >>>         signal.sepfir2d(ck, [1], derfilt)

   Alternatively we could have done::

       laplacian = array([[0,1,0],[1,-4,1],[0,1,0]],float32)
       deriv2 = signal.convolve2d(ck,laplacian,mode='same',boundary='symm')

   >>> plt.figure()
   >>> plt.imshow(image)
   >>> plt.gray()
   >>> plt.title('Original image')
   >>> plt.show()

   >>> plt.figure()
   >>> plt.imshow(deriv)
   >>> plt.gray()
   >>> plt.title('Output of spline edge filter')
   >>> plt.show()

..   :caption: Example of using smoothing splines to filter images.


Filtering
---------

Filtering is a generic name for any system that modifies an input
signal in some way. In SciPy a signal can be thought of as a Numpy
array. There are different kinds of filters for different kinds of
operations. There are two broad kinds of filtering operations: linear
and non-linear. Linear filters can always be reduced to multiplication
of the flattened Numpy array by an appropriate matrix resulting in
another flattened Numpy array. Of course, this is not usually the best
way to compute the filter as the matrices and vectors involved may be
huge. For example filtering a :math:`512 \times 512` image with this
method would require multiplication of a :math:`512^2 \times 512^2`
matrix with a :math:`512^2` vector. Just trying to store the
:math:`512^2 \times 512^2` matrix using a standard Numpy array would
require :math:`68,719,476,736` elements. At 4 bytes per element this
would require :math:`256\textrm{GB}` of memory. In most applications
most of the elements of this matrix are zero and a different method
for computing the output of the filter is employed.


Convolution/Correlation
^^^^^^^^^^^^^^^^^^^^^^^

Many linear filters also have the property of shift-invariance. This
means that the filtering operation is the same at different locations
in the signal and it implies that the filtering matrix can be
constructed from knowledge of one row (or column) of the matrix alone.
In this case, the matrix multiplication can be accomplished using
Fourier transforms.

Let :math:`x\left[n\right]` define a one-dimensional signal indexed by the
integer :math:`n.` Full convolution of two one-dimensional signals can be
expressed as

.. math::
   :nowrap:

    \[ y\left[n\right]=\sum_{k=-\infty}^{\infty}x\left[k\right]h\left[n-k\right].\]

This equation can only be implemented directly if we limit the
sequences to finite support sequences that can be stored in a
computer, choose :math:`n=0` to be the starting point of both
sequences, let :math:`K+1` be that value for which
:math:`y\left[n\right]=0` for all :math:`n>K+1` and :math:`M+1` be
that value for which :math:`x\left[n\right]=0` for all :math:`n>M+1` ,
then the discrete convolution expression is

.. math::
   :nowrap:

    \[ y\left[n\right]=\sum_{k=\max\left(n-M,0\right)}^{\min\left(n,K\right)}x\left[k\right]h\left[n-k\right].\]

For convenience assume :math:`K\geq M.` Then, more explicitly the output of
this operation is

.. math::
   :nowrap:

    \begin{eqnarray*} y\left[0\right] & = & x\left[0\right]h\left[0\right]\\ y\left[1\right] & = & x\left[0\right]h\left[1\right]+x\left[1\right]h\left[0\right]\\ y\left[2\right] & = & x\left[0\right]h\left[2\right]+x\left[1\right]h\left[1\right]+x\left[2\right]h\left[0\right]\\ \vdots & \vdots & \vdots\\ y\left[M\right] & = & x\left[0\right]h\left[M\right]+x\left[1\right]h\left[M-1\right]+\cdots+x\left[M\right]h\left[0\right]\\ y\left[M+1\right] & = & x\left[1\right]h\left[M\right]+x\left[2\right]h\left[M-1\right]+\cdots+x\left[M+1\right]h\left[0\right]\\ \vdots & \vdots & \vdots\\ y\left[K\right] & = & x\left[K-M\right]h\left[M\right]+\cdots+x\left[K\right]h\left[0\right]\\ y\left[K+1\right] & = & x\left[K+1-M\right]h\left[M\right]+\cdots+x\left[K\right]h\left[1\right]\\ \vdots & \vdots & \vdots\\ y\left[K+M-1\right] & = & x\left[K-1\right]h\left[M\right]+x\left[K\right]h\left[M-1\right]\\ y\left[K+M\right] & = & x\left[K\right]h\left[M\right].\end{eqnarray*}

Thus, the full discrete convolution of two finite sequences of lengths
:math:`K+1` and :math:`M+1` respectively results in a finite sequence of length
:math:`K+M+1=\left(K+1\right)+\left(M+1\right)-1.`

One dimensional convolution is implemented in SciPy with the function
``signal.convolve`` . This function takes as inputs the signals
:math:`x,` :math:`h` , and an optional flag and returns the signal
:math:`y.` The optional flag allows for specification of which part of
the output signal to return. The default value of 'full' returns the
entire signal. If the flag has a value of 'same' then only the middle
:math:`K` values are returned starting at :math:`y\left[\left\lfloor
\frac{M-1}{2}\right\rfloor \right]` so that the output has the same
length as the largest input. If the flag has a value of 'valid' then
only the middle :math:`K-M+1=\left(K+1\right)-\left(M+1\right)+1`
output values are returned where :math:`z` depends on all of the
values of the smallest input from :math:`h\left[0\right]` to
:math:`h\left[M\right].` In other words only the values
:math:`y\left[M\right]` to :math:`y\left[K\right]` inclusive are
returned.

This same function ``signal.convolve`` can actually take :math:`N`
-dimensional arrays as inputs and will return the :math:`N`
-dimensional convolution of the two arrays. The same input flags are
available for that case as well.

Correlation is very similar to convolution except for the minus sign
becomes a plus sign. Thus

.. math::
   :nowrap:

    \[ w\left[n\right]=\sum_{k=-\infty}^{\infty}y\left[k\right]x\left[n+k\right]\]

is the (cross) correlation of the signals :math:`y` and :math:`x.` For
finite-length signals with :math:`y\left[n\right]=0` outside of the range
:math:`\left[0,K\right]` and :math:`x\left[n\right]=0` outside of the range
:math:`\left[0,M\right],` the summation can simplify to

.. math::
   :nowrap:

    \[ w\left[n\right]=\sum_{k=\max\left(0,-n\right)}^{\min\left(K,M-n\right)}y\left[k\right]x\left[n+k\right].\]

Assuming again that :math:`K\geq M` this is

.. math::
   :nowrap:

    \begin{eqnarray*} w\left[-K\right] & = & y\left[K\right]x\left[0\right]\\ w\left[-K+1\right] & = & y\left[K-1\right]x\left[0\right]+y\left[K\right]x\left[1\right]\\ \vdots & \vdots & \vdots\\ w\left[M-K\right] & = & y\left[K-M\right]x\left[0\right]+y\left[K-M+1\right]x\left[1\right]+\cdots+y\left[K\right]x\left[M\right]\\ w\left[M-K+1\right] & = & y\left[K-M-1\right]x\left[0\right]+\cdots+y\left[K-1\right]x\left[M\right]\\ \vdots & \vdots & \vdots\\ w\left[-1\right] & = & y\left[1\right]x\left[0\right]+y\left[2\right]x\left[1\right]+\cdots+y\left[M+1\right]x\left[M\right]\\ w\left[0\right] & = & y\left[0\right]x\left[0\right]+y\left[1\right]x\left[1\right]+\cdots+y\left[M\right]x\left[M\right]\\ w\left[1\right] & = & y\left[0\right]x\left[1\right]+y\left[1\right]x\left[2\right]+\cdots+y\left[M-1\right]x\left[M\right]\\ w\left[2\right] & = & y\left[0\right]x\left[2\right]+y\left[1\right]x\left[3\right]+\cdots+y\left[M-2\right]x\left[M\right]\\ \vdots & \vdots & \vdots\\ w\left[M-1\right] & = & y\left[0\right]x\left[M-1\right]+y\left[1\right]x\left[M\right]\\ w\left[M\right] & = & y\left[0\right]x\left[M\right].\end{eqnarray*}



The SciPy function ``signal.correlate`` implements this
operation. Equivalent flags are available for this operation to return
the full :math:`K+M+1` length sequence ('full') or a sequence with the
same size as the largest sequence starting at
:math:`w\left[-K+\left\lfloor \frac{M-1}{2}\right\rfloor \right]`
('same') or a sequence where the values depend on all the values of
the smallest sequence ('valid'). This final option returns the
:math:`K-M+1` values :math:`w\left[M-K\right]` to
:math:`w\left[0\right]` inclusive.

The function :obj:`signal.correlate` can also take arbitrary :math:`N`
-dimensional arrays as input and return the :math:`N` -dimensional
convolution of the two arrays on output.

When :math:`N=2,` :obj:`signal.correlate` and/or
:obj:`signal.convolve` can be used to construct arbitrary image
filters to perform actions such as blurring, enhancing, and
edge-detection for an image.

Convolution is mainly used for filtering when one of the signals is
much smaller than the other ( :math:`K\gg M` ), otherwise linear
filtering is more easily accomplished in the frequency domain (see
Fourier Transforms).


Difference-equation filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A general class of linear one-dimensional filters (that includes
convolution filters) are filters described by the difference equation

.. math::
   :nowrap:

    \[ \sum_{k=0}^{N}a_{k}y\left[n-k\right]=\sum_{k=0}^{M}b_{k}x\left[n-k\right]\]

where :math:`x\left[n\right]` is the input sequence and
:math:`y\left[n\right]` is the output sequence. If we assume initial
rest so that :math:`y\left[n\right]=0` for :math:`n<0` , then this
kind of filter can be implemented using convolution.  However, the
convolution filter sequence :math:`h\left[n\right]` could be infinite
if :math:`a_{k}\neq0` for :math:`k\geq1.` In addition, this general
class of linear filter allows initial conditions to be placed on
:math:`y\left[n\right]` for :math:`n<0` resulting in a filter that
cannot be expressed using convolution.

The difference equation filter can be thought of as finding
:math:`y\left[n\right]` recursively in terms of it's previous values

.. math::
   :nowrap:

    \[ a_{0}y\left[n\right]=-a_{1}y\left[n-1\right]-\cdots-a_{N}y\left[n-N\right]+\cdots+b_{0}x\left[n\right]+\cdots+b_{M}x\left[n-M\right].\]

Often :math:`a_{0}=1` is chosen for normalization. The implementation
in SciPy of this general difference equation filter is a little more
complicated then would be implied by the previous equation. It is
implemented so that only one signal needs to be delayed. The actual
implementation equations are (assuming :math:`a_{0}=1` ).

.. math::
   :nowrap:

    \begin{eqnarray*} y\left[n\right] & = & b_{0}x\left[n\right]+z_{0}\left[n-1\right]\\ z_{0}\left[n\right] & = & b_{1}x\left[n\right]+z_{1}\left[n-1\right]-a_{1}y\left[n\right]\\ z_{1}\left[n\right] & = & b_{2}x\left[n\right]+z_{2}\left[n-1\right]-a_{2}y\left[n\right]\\ \vdots & \vdots & \vdots\\ z_{K-2}\left[n\right] & = & b_{K-1}x\left[n\right]+z_{K-1}\left[n-1\right]-a_{K-1}y\left[n\right]\\ z_{K-1}\left[n\right] & = & b_{K}x\left[n\right]-a_{K}y\left[n\right],\end{eqnarray*}

where :math:`K=\max\left(N,M\right).` Note that :math:`b_{K}=0` if
:math:`K>M` and :math:`a_{K}=0` if :math:`K>N.` In this way, the
output at time :math:`n` depends only on the input at time :math:`n`
and the value of :math:`z_{0}` at the previous time. This can always
be calculated as long as the :math:`K` values
:math:`z_{0}\left[n-1\right]\ldots z_{K-1}\left[n-1\right]` are
computed and stored at each time step.

The difference-equation filter is called using the command
:obj:`signal.lfilter` in SciPy. This command takes as inputs the
vector :math:`b,` the vector, :math:`a,` a signal :math:`x` and
returns the vector :math:`y` (the same length as :math:`x` ) computed
using the equation given above. If :math:`x` is :math:`N`
-dimensional, then the filter is computed along the axis provided. If,
desired, initial conditions providing the values of
:math:`z_{0}\left[-1\right]` to :math:`z_{K-1}\left[-1\right]` can be
provided or else it will be assumed that they are all zero. If initial
conditions are provided, then the final conditions on the intermediate
variables are also returned. These could be used, for example, to
restart the calculation in the same state.

Sometimes it is more convenient to express the initial conditions in
terms of the signals :math:`x\left[n\right]` and
:math:`y\left[n\right].` In other words, perhaps you have the values
of :math:`x\left[-M\right]` to :math:`x\left[-1\right]` and the values
of :math:`y\left[-N\right]` to :math:`y\left[-1\right]` and would like
to determine what values of :math:`z_{m}\left[-1\right]` should be
delivered as initial conditions to the difference-equation filter. It
is not difficult to show that for :math:`0\leq m<K,`

.. math::
   :nowrap:

    \[ z_{m}\left[n\right]=\sum_{p=0}^{K-m-1}\left(b_{m+p+1}x\left[n-p\right]-a_{m+p+1}y\left[n-p\right]\right).\]

Using this formula we can find the intial condition vector
:math:`z_{0}\left[-1\right]` to :math:`z_{K-1}\left[-1\right]` given initial
conditions on :math:`y` (and :math:`x` ). The command :obj:`signal.lfiltic`
performs this function.


Other filters
^^^^^^^^^^^^^

The signal processing package provides many more filters as well.


Median Filter
"""""""""""""

A median filter is commonly applied when noise is markedly non-Gaussian
or when it is desired to preserve edges. The median filter
works by sorting all of the array pixel values in a rectangular region
surrounding the point of interest. The sample median of this list of
neighborhood pixel values is used as the value for the output array.
The sample median is the middle array value in a sorted list of
neighborhood values. If there are an even number of elements in the
neighborhood, then the average of the middle two values is used as the
median. A general purpose median filter that works on N-dimensional
arrays is :obj:`signal.medfilt` . A specialized version that works
only for two-dimensional arrays is available as
:obj:`signal.medfilt2d` .


Order Filter
""""""""""""

A median filter is a specific example of a more general class of
filters called order filters. To compute the output at a particular
pixel, all order filters use the array values in a region surrounding
that pixel. These array values are sorted and then one of them is
selected as the output value. For the median filter, the sample median
of the list of array values is used as the output. A general order
filter allows the user to select which of the sorted values will be
used as the output. So, for example one could choose to pick the
maximum in the list or the minimum. The order filter takes an
additional argument besides the input array and the region mask that
specifies which of the elements in the sorted list of neighbor array
values should be used as the output. The command to perform an order
filter is :obj:`signal.order_filter` .


Wiener filter
"""""""""""""

The Wiener filter is a simple deblurring filter for denoising images.
This is not the Wiener filter commonly described in image
reconstruction problems but instead it is a simple, local-mean filter.
Let :math:`x` be the input signal, then the output is



.. math::
   :nowrap:

    \[ y=\left\{ \begin{array}{cc} \frac{\sigma^{2}}{\sigma_{x}^{2}}m_{x}+\left(1-\frac{\sigma^{2}}{\sigma_{x}^{2}}\right)x & \sigma_{x}^{2}\geq\sigma^{2},\\ m_{x} & \sigma_{x}^{2}<\sigma^{2},\end{array}\right.\]

where :math:`m_{x}` is the local estimate of the mean and
:math:`\sigma_{x}^{2}` is the local estimate of the variance. The
window for these estimates is an optional input parameter (default is
:math:`3\times3` ). The parameter :math:`\sigma^{2}` is a threshold
noise parameter. If :math:`\sigma` is not given then it is estimated
as the average of the local variances.


Hilbert filter
""""""""""""""

The Hilbert transform constructs the complex-valued analytic signal
from a real signal. For example if :math:`x=\cos\omega n` then
:math:`y=\textrm{hilbert}\left(x\right)` would return (except near the
edges) :math:`y=\exp\left(j\omega n\right).` In the frequency domain,
the hilbert transform performs

.. math::
   :nowrap:

    \[ Y=X\cdot H\]

where :math:`H` is 2 for positive frequencies, :math:`0` for negative
frequencies and :math:`1` for zero-frequencies.


Least-Squares Spectral Analysis (:mod:`spectral`)
-------------------------------------------------

Least-squares spectral analysis (LSSA) is a method of estimating a frequency
spectrum, based on a least squares fit of sinusoids to data samples, similar to
Fourier analysis. Fourier analysis, the most used spectral method in science,
generally boosts long-periodic noise in long gapped records; LSSA mitigates
such problems.


Lomb-Scargle Periodograms (:func:`spectral.lombscargle`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Lomb-Scargle method performs spectral analysis on unevenly sampled data and
is known to be a powerful way to find, and test the significance of, weak
periodic signals.

For a time series comprising :math:`N_{t}` measurements
:math:`X_{j}\equiv X(t_{j})` sampled at times :math:`t_{j}` where
:math:`(j = 1, \ldots, N_{t})`, assumed to have been scaled and shifted
such that its mean is zero and its variance is unity, the normalized
Lomb-Scargle periodogram at frequency :math:`f` is

.. math::

    P_{n}(f) \frac{1}{2}\left\{\frac{\left[\sum_{j}^{N_{t}}X_{j}\cos\omega(t_{j}-\tau)\right]^{2}}{\sum_{j}^{N_{t}}\cos^{2}\omega(t_{j}-\tau)}+\frac{\left[\sum_{j}^{N_{t}}X_{j}\sin\omega(t_{j}-\tau)\right]^{2}}{\sum_{j}^{N_{t}}\sin^{2}\omega(t_{j}-\tau)}\right\}.

Here, :math:`\omega \equiv 2\pi f` is the angular frequency.
The frequency dependent time offset :math:`\tau` is given by

.. math::

    \tan 2\omega\tau = \frac{\sum_{j}^{N_{t}}\sin 2\omega t_{j}}{\sum_{j}^{N_{t}}\cos 2\omega t_{j}}.

The :func:`~scipy.signal.spectral.lombscargle` function
calculates the periodogram using a slightly
modified algorithm due to Townsend [3]_ which allows the
periodogram to be calculated using only a single pass through
the input arrays for each frequency.

The equation is refactored as:

.. math::

    P_{n}(f) = \frac{1}{2}\left[\frac{(c_{\tau}XC + s_{\tau}XS)^{2}}{c_{\tau}^{2}CC + 2c_{\tau}s_{\tau}CS + s_{\tau}^{2}SS} + \frac{(c_{\tau}XS - s_{\tau}XC)^{2}}{c_{\tau}^{2}SS - 2c_{\tau}s_{\tau}CS + s_{\tau}^{2}CC}\right]

and

.. math::

    \tan 2\omega\tau = \frac{2CS}{CC-SS}.

Here,

.. math::

    c_{\tau} = \cos\omega\tau,\qquad s_{\tau} = \sin\omega\tau

while the sums are

.. math::

    XC &= \sum_{j}^{N_{t}} X_{j}\cos\omega t_{j}\\
    XS &= \sum_{j}^{N_{t}} X_{j}\sin\omega t_{j}\\
    CC &= \sum_{j}^{N_{t}} \cos^{2}\omega t_{j}\\
    SS &= \sum_{j}^{N_{t}} \sin^{2}\omega t_{j}\\
    CS &= \sum_{j}^{N_{t}} \cos\omega t_{j}\sin\omega t_{j}.

This requires :math:`N_{f}(2N_{t}+3)` trigonometric function
evaluations giving a factor of :math:`\sim 2` speed increase over the
straightforward implementation.


.. XXX: TODO
..
.. Detrend
.. """""""
..
.. Filter design
.. -------------
..
..
.. Finite-impulse response design
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Inifinite-impulse response design
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Analog filter frequency response
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Digital filter frequency response
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Linear Time-Invariant Systems
.. -----------------------------
..
..
.. LTI Object
.. ^^^^^^^^^^
..
..
.. Continuous-Time Simulation
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Step response
.. ^^^^^^^^^^^^^
..
..
.. Impulse response
.. ^^^^^^^^^^^^^^^^
..
..
.. Input/Output
.. ============
..
..
.. Binary
.. ------
..
..
.. Arbitrary binary input and output (fopen)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Read and write Matlab .mat files
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Saving workspace
.. ^^^^^^^^^^^^^^^^
..
..
.. Text-file
.. ---------
..
..
.. Read text-files (read_array)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Write a text-file (write_array)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..
..
.. Fourier Transforms
.. ==================
..
..
.. One-dimensional
.. ---------------
..
..
.. Two-dimensional
.. ---------------
..
..
.. N-dimensional
.. -------------
..
..
.. Shifting
.. --------
..
..
.. Sample frequencies
.. ------------------
..
..
.. Hilbert transform
.. -----------------
..
..
.. Tilbert transform
.. -----------------


.. rubric:: References

Some further reading and related software:

.. [1] N.R. Lomb "Least-squares frequency analysis of unequally spaced
       data", Astrophysics and Space Science, vol 39, pp. 447-462, 1976

.. [2] J.D. Scargle "Studies in astronomical time series analysis. II - 
       Statistical aspects of spectral analysis of unevenly spaced data",
       The Astrophysical Journal, vol 263, pp. 835-853, 1982

.. [3] R.H.D. Townsend, "Fast calculation of the Lomb-Scargle
       periodogram using graphics processing units.", The Astrophysical
       Journal Supplement Series, vol 191, pp. 247-253, 2010

