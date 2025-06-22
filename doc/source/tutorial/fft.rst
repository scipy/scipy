Fourier Transforms (:mod:`scipy.fft`)
=========================================

.. sectionauthor:: SciPy Developers

.. currentmodule:: scipy.fft

.. contents::


Fourier analysis refers to techniques for decomposing functions into periodic
components as well as synthesising functions from such periodic components. The
"Discrete Fourier Transform" (DFT) calculates discrete frequency components from a
sequence of equidistant samples representing an input function. The DFT has become a
mainstay of numerical computing in part because it can be efficiently implemented with
a complexity of :math:`O(n \log n)` by an algorithm called the "Fast Fourier Transform"
(FFT). This algorithm was known to Gauss (1805) and was brought to light in its current
form by Cooley and Tukey [CT65]_. Press et al. [NR07]_ provide an accessible
introduction to Fourier analysis and its applications.

.. _tutorial_FFT:

Fast Fourier transforms
-----------------------

1-D DFT
_______

The DFT of a complex-valued sequence :math:`x[n]` of length :math:`N,` is defined as

.. math::

    y[k] = \sum_{n=0}^{N-1} e^{-j 2 \pi k n / N} x[n] \, ,

where :math:`y[k]` is also a complex-valued sequence of length :math:`N`. Due to the
periodicity of the complex exponential, the DFT is :math:`N`-periodic, i.e.,
:math:`y[k+lN] = y[k]`, :math:`l\in\mathbb{Z}`. The inverse transform is given by

.. math::

    x[n] = \frac{1}{N} \sum_{k=0}^{N-1} e^{j 2 \pi k n / N} y[k] \, .

Note that the :math:`N`-periodicity holds also for inverse DFT values, i.e.,
:math:`x[n+lN] = x[n]`, :math:`l\in\mathbb{Z}`.

Since the DFT is a linear transform in :math:`\mathbb{C}^N`, it can be expressed as a
symmetric :math:`N\times N` transformation matrix :math:`\mathbf{F}` with elements
:math:`F[k,n] = e^{-j 2 \pi k n / N}`. The transformation matrix of the inverse DFT
:math:`\mathbf{F}^{-1} = \overline{\mathbf{F}^T}/N` is the conjugate transpose of
:math:`\mathbf{F}` divided by :math:`N`.

It can be shown that

.. math::

    y[0] = \sum_{n=0}^{N-1} x[n] \ ,\quad
    x[0] = \frac{1}{N} \sum_{k=0}^{N-1} y[k] \ ,\quad
    \sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |y[k]|^2

The DFT and its inverse can be calculated by utilizing the :func:`fft`
and :func:`ifft` functions, as illustrated in the following example:

>>> import numpy as np
>>> from scipy.fft import fft, ifft
...
>>> x = np.array([2. , 1.5, 0.5, 0. , 0.5, 1.5])
>>> y = fft(x)
>>> y
array([6.-0.j, 3.+0.j, 0.+0.j, 0.-0.j, 0.-0.j, 3.-0.j])
>>> ifft(y)  # equals x, but is complex-valued
array([2. +0.j, 1.5+0.j, 0.5+0.j, 0. +0.j, 0.5-0.j, 1.5+0.j])
>>> # the following are equal for any complex x:
>>> np.isclose(ifft(x), fft(x.conj()).conj()/len(x))
array([ True,  True,  True,  True,  True,  True])
>>> y[0], np.sum(x)
(np.complex128(6-0j), np.float64(6.0))
>>> x[0], np.mean(y)
(np.float64(2.0), np.complex128(2+0j))
>>> np.sum(np.abs(x)**2), np.mean(np.abs(y)**2)
(np.float64(9.0), np.float64(9.0))

The following example illustrates the connection between the Fourier series

.. math::

    x(t) = \sum_{k\in\mathbb{Z}} c_k\, e^{j 2 \pi k t}

and the DFT. :math:`x(t)` may be interpreted as a two-dimensional parametric closed
curve :math:`x: t\in[0, 1] \mapsto \mathbb{C}`, where the real and the imaginary part
of :math:`x(t)` represent represent the two dimensions. This curve is sampled at
equistant intervals, i.e., :math:`x[n] := x(n/N)` and the Fourier coefficents
:math:`c_k = y[k] / N` can be efficiently calculated using a DFT. Here, a curve,
representing the greek letter ω, is generated from 17 non-zero DFT values. The
:func:`ifft`
function is used to calculate the 31 sample values :math:`x[n]` and a Fourier series is
used to interpolate between those samples:

.. plot:: tutorial/examples/fft_ParemetricClosedCurve.py
    :alt: This code generates from DFT values a plot of a parametric closed curve
          representing the greek letter omega.

The samples are depicted as dots and the value of the parameter :math:`t` is encoded in
the color of the curve. Note that the `~scipy.signal.resample` function, which uses the
:func:`fft`/:func:`ifft` functions internally, can be used as a more efficient
alternative for calculating the Fourier expansion for each interpolation value.


In the usual case in which :math:`x[n]` is real, the subsequence
:math:`y[1],\ldots,y[N-1]` is conjugate-symmetric, i.e.,
:math:`y[N-k] = y^*[k],` for :math:`k=1,\ldots,N-1.`
Therefore :math:`y[N/2]` is real for even :math:`N,`
while :math:`y[0]` is real for either even or odd :math:`N.`
The subsequence :math:`y[0],\ldots,y[N//2]` contains :math:`N`
independent real values, and it is sufficient to reconstruct :math:`x[n].`
Typically only this subsequence is plotted for real :math:`x[n].`

In many applications, the sequence :math:`x[n]` represents samples taken
from a real function of time :math:`x(t),` at a fixed time interval :math:`T_s,`
i.e., :math:`x[n] = x(n\,T_s).` The parameter :math:`T_s` is called
"sampling time", or "time resolution", while its inverse :math:`F_s=1/T_s`
is called "sampling frequency", or "sampling rate".
The time duration spanned by :math:`N` samples is :math:`T_p=N\,T_s`
and is called "sampling period", while its inverse
:math:`F_p=1/T_p=F_s/N` is called "bin width", or "frequency resolution".
In analogy to :math:`x[n],` the sequence :math:`y[k]` is considered to be
the samples of a complex function of frequency :math:`y(f),` at a fixed
frequency interval :math:`F_p,` i.e., :math:`y[k] = y(k\,F_p).`
The frequencies :math:`k\,F_p` are called "bin frequencies".

The following example models the samples taken from a real signal
containing two sinusoidal components and a bias (such that its
spectrum can be easily predicted). It calculates the DFT of the signal
using :func:`fft`, and then it plots its amplitude vs. frequency.
The frequency points are caclulated using the convenience function
:func:`rfftfreq`, which takes as arguments the parameters :math:`N`
and :math:`T_s` and returns the multiples of the bin width :math:`F_p`
in the range :math:`[0,F_s/2].` These are the frequency points of
interest for the FFT of a real signal.

.. plot::
    :alt: "This code generates a plot of the amplitude of the transform vs. frequency, shown as a blue line. Its Y value is zero all the way across the X axis, with the exception of three peaks. The first peak is at zero frequency, with expected amplitude 0.4. The second peak is at 50 Hz with expected amplitude 0.8. The third peak is at 125 Hz with expected amplitude 0.6. The expected amplitudes at each frequency are shown with short horizontal lines. More details around the peaks at 50 Hz and 125 Hz are shown in inset axes."

    >>> from scipy.fft import fft, rfftfreq
    >>> import numpy as np
    >>> N = 400        # number of samples
    >>> Fs = 800.0     # sampling frequency in Hz
    >>> Ts = 1.0 / Fs  # sampling time (step) in s
    >>> Tp = N / Fs    # sampling period (total) in s
    >>> Fp = Fs / N    # bin width (frequency resolution) in Hz
    >>> Fx = np.array([0.0, 50.0, 125.0])            # signal frequencies in Hz
    >>> Ax = np.array([0.2, 0.8, 0.6])               # signal amplitudes
    >>> t = np.linspace(0.0, Tp, N, endpoint=False)  # time in s
    >>> f = rfftfreq(N, Ts)                          # frequency in Hz
    >>> x = np.sum(Ax[:,np.newaxis]*np.cos(2.0*np.pi*np.outer(Fx,t)), axis=0)
    >>> y = fft(x)
    >>> ya = 2.0/N * np.abs(y[:len(f)])
    >>>
    >>> # plot ya vs. f in the range [0, Fs/2]
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(7,4))
    >>> fig.tight_layout()
    >>> ax.plot(f, ya)
    >>>
    >>> # decorations
    >>> for i in range(len(Fx)):
    ...     # expected amplitude at frequency Fx[i]
    ...     a = Ax[i] if i > 0 else 2*Ax[i]
    ...     ax.plot(Fx[i]+np.array([[-6.0, 6.0], [-3.0, 3.0]]), a*np.r_[1.0, 1.0], 'g')
    >>> ax.set_title(f'N={N}, Fp={Fp} Hz, Fs={Fs} Hz')
    >>> ax.set_xlabel('f [Hz]')
    >>> ax.set_ylabel('2/N |y|')
    >>> ax.grid()
    >>> for i in (1, 2):
    ...     axi = fig.add_axes([0.35+0.22*i, 0.4, 0.13, 0.4])
    ...     # frequency indices around Fx[i]
    ...     k = np.r_[int(Fx[i]/Fp+0.3)-4:int(Fx[i]/Fp+0.8)+5]
    ...     axi.scatter(f[k], ya[k], marker='o')
    ...     axi.plot(np.outer([1, 1], f[k]), np.outer([0, 1], ya[k]), 'b')
    ...     axi.set_xlim(f[k[0]]-Fp/2, f[k[-1]]+Fp/2)
    ...     axi.set_ylim(-0.1, 0.9)
    ...     axi.grid()
    >>> plt.show()


The plotted spectrum is scaled such that the value at a peak gives
the amplitude of the sinusoidal component at that frequency.
With this scaling, the peak at zero frequency is equal to the mean value
of the signal multiplied by 2.
The first sinusoidal component has frequency equal to an exact multiple
of :math:`F_p,` i.e., equal to one of the bin frequencies, therefore it
appears as an impulse with the expected amplitude.
The whole energy (the square of the Euclidean norm) of this component
is concentated in a single bin. However, the second sinusoidal component
has a frequency between two bin frequencies and its energy is spread in
multiple bins (theoretically in all bins).
This effect is called "spectral leakage" (see Appendix).
Taking this effect into account, the expected amplitude at each of
the two main peaks is :math:`2/\pi\cdot A_x[2],` i.e., 81% of the
component energy is contained in the two main peaks and 19% is spread
over the whole frequency range.

In the previous example we calculated the DFT of a sequence :math:`x`
using :func:`fft`, and then we scaled its amplitude (absolute value)
by :math:`2/N,` as in the following code:

>>> y = fft(x)
>>> ya = 2.0/N * np.abs(y[:len(f)])

Since :math:`x` is real and the phase (angle) of :math:`y` is of no interest,
we can obtain the same result by calculating the inverse DFT of :math:`x`
using :func:`ifft`, and then scaling its amplitude by :math:`2.`

>>> yi = ifft(x)
>>> ya = 2.0 * np.abs(yi[:len(f)])

For a real sequence :math:`x` of length :math:`N,` only the first
:math:`N//2+1` values of its DFT are used, since the remaining
:math:`(N-1)//2` values are redundant.
The function :func:`rfft` calculates the FFT of a real sequence
and returns only the first :math:`N//2+1` values, which correspond to
the :math:`N//2+1` frequency values returned by :func:`rfftfreq`.

>>> from scipy.fft import rfft
>>> yr = rfft(x)
>>> ya = 2.0/N * np.abs(yr)

The function :func:`irfft` is the inverse of :func:`rfft`, i.e.,
it calculates the real sequence :math:`x[0],\ldots,x[N-1]` from
a subsequence :math:`y[0],\ldots,y[N//2]` of its DFT.
The parameter :math:`N` must be specified as an additional argument ``n``,
because its parity cannot be determined from the length :math:`N//2+1`
of the input argument (if ``n`` is not specified, an even :math:`N` is assumed).

>>> from scipy.fft import fft, rfft, irfft
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5, 1.0])
>>> fft(x)
array([ 5.5 +0.j        ,  2.25-0.4330127j , -2.75-1.29903811j,
        1.5 +0.j        , -2.75+1.29903811j,  2.25+0.4330127j ])
>>> yr = rfft(x)
>>> yr
array([ 5.5 +0.j        ,  2.25-0.4330127j , -2.75-1.29903811j,
        1.5 +0.j        ])
>>> irfft(yr)
array([ 1. ,  2. ,  1. , -1. ,  1.5,  1. ])
>>>
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])
>>> fft(x)
array([ 4.5       +0.j        ,  2.08155948-1.65109876j,
       -1.83155948+1.60822041j, -1.83155948-1.60822041j,
        2.08155948+1.65109876j])
>>> yr = rfft(x)
>>> yr
array([ 4.5       +0.j        ,  2.08155948-1.65109876j,
        -1.83155948+1.60822041j])
>>> # n must be specified if it is odd
>>> irfft(yr, n=len(x))
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>> # otherwise the result is not equal to the original signal
>>> irfft(yr)
array([ 1.70788987,  2.40843925, -0.37366961,  0.75734049])


The DFT of a sequence :math:`x` of length :math:`N,` can be extended to
an infinite sequence :math:`y[k]` by using the DFT formula for any integer
:math:`k.` With this extension, the infinite sequence :math:`y[k]` is a
periodic function of :math:`k` with period :math:`N,` and the index :math:`k`
is considered to be wrapped with period :math:`N,` i.e., the values
:math:`k+mN` for different integers :math:`m` are indistinguishable.
Similarly, the bin frequency in the DFT of a sampled signal is considered
to be wrapped with period equal to the sampling frequency :math:`F_s.`

If the sequence :math:`x` is not real, then the subsequence
:math:`y[1],\ldots,y[N-1]` of its DFT is not conjugate-symmetric.
The corresponding amplitudes may be symmetric (by chance or
by construction), but usually they are not.
In this case it is common to plot the spectrum of a signal with the
zero frequency at the center of the horizontal axis, with the
understanding that the bin frequency is wrapped with period :math:`F_s.`
For an even :math:`N,` the bin frequency :math:`F_s/2` can be shown
either as positive or negative (but usually not both).
This kind of plot is called "two-sided", while the plot of the spectrum
of a real signal in the range :math:`[0,F_s/2]` is called "one-sided".

To simplify working with a two-sided spectrum, scipy provides two helper
functions. The function :func:`fftfreq` takes as arguments the parameters
:math:`N` and :math:`T_s,` and returns the sequence of bin frequencies
of length :math:`N,` in which all bin frequencies :math:`f \ge F_s/2`
are wrapped to :math:`f-F_s.` This makes negative the last :math:`N//2`
bin frequencies in the sequence. The companion function :func:`fftshift`
takes as argument a sequence of length :math:`N,` and returns the same
sequence rotated to the right by :math:`N//2.` When it is applied to
the output of :func:`fftfreq`, the result is in arithmetic order,
as needed for the plot of a two-sided spectrum.

>>> from scipy.fft import rfftfreq, fftfreq, fftshift
>>> N = 8
>>> rfftfreq(N, 1.0/N)
array([0., 1., 2., 3., 4.])
>>> f = fftfreq(N, 1.0/N)
>>> f
array([ 0., 1., 2., 3., -4., -3., -2., -1.])
>>> fftshift(f)
array([-4., -3., -2., -1.,  0.,  1.,  2.,  3.])
>>> x = np.arange(N)
>>> fftshift(x)
array([4, 5, 6, 7, 0, 1, 2, 3])

Notice that for an even :math:`N,` as in the example above,
the bin frequency :math:`F_s/2` is shown as positive in a one-sided
spectrum created with :func:`rfftfreq`, and as negative in a two-sided
spectrum created with :func:`fftfreq` and :func:`fftshift`.

The example below plots the spectrum of a complex signal containing
two complex exponential components.

.. plot::
    :alt: "This code generates an X-Y plot with amplitude on the Y axis vs. frequency on the X axis. The trace is zero-valued across the plot except for two sharp peaks at -80 and 50 Hz. The 50 Hz peak on the right is twice as tall."

    >>> from scipy.fft import fft, fftfreq, fftshift
    >>> import numpy as np
    >>> N = 400        # number of samples
    >>> Fs = 800.0     # sampling frequency in Hz
    >>> Ts = 1.0 / Fs  # sampling time (step) in s
    >>> Tp = N / Fs    # sampling period (total) in s
    >>> Fp = Fs / N    # bin width (frequency resolution) in Hz
    >>> t = np.linspace(0.0, Tp, N, endpoint=False)  # time in s
    >>> f = fftfreq(N, Ts)                           # frequency in Hz
    >>> x = 0.8*np.exp(1.j * 2.0*np.pi*50.0*t) + 0.4*np.exp(-1.j * 2.0*np.pi*80.0*t)
    >>> y = fft(x)
    >>> ya = 1.0/N * np.abs(y)
    >>> fc = fftshift(f)
    >>> yc = fftshift(ya)
    >>>
    >>> # plot yc vs. fc in the range [-Fs/2, Fs/2)
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(6,4))
    >>> fig.tight_layout()
    >>> ax.plot(fc, yc)
    >>>
    >>> # decorations
    >>> ax.plot( 50.0+np.array([[-12.0, 12.0], [-8.0, 8.0]]), 0.8*np.r_[1.0, 1.0], 'g')
    >>> ax.plot(-80.0+np.array([[-12.0, 12.0], [-8.0, 8.0]]), 0.4*np.r_[1.0, 1.0], 'g')
    >>> ax.set_ylim(-0.1, 0.9)
    >>> ax.set_title(f'N={N}, Fp={Fp} Hz, Fs={Fs} Hz')
    >>> ax.set_xlabel('f [Hz]')
    >>> ax.set_ylabel('|y|')
    >>> ax.grid()
    >>> plt.show()

Notice that the spectrum is asymmetric with respect to zero frequency.
If it had been plotted in the range :math:`[0,F_s),`
the peak at -80 Hz would appear at 720 Hz.
The amplitude is scaled such that the value at a peak gives
the amplitude of the complex exponential component at that frequency.
No spectral leakage is visible, since both signal frequencies are
exact multiples of the bin width.


2-D and N-D DFT
_______________

The functions :func:`fft2` and :func:`ifft2` provide 2-D FFT and
IFFT, respectively. Similarly, :func:`fftn` and :func:`ifftn` provide
N-D FFT, and IFFT, respectively.

For real-input signals, similarly to :func:`rfft`, we have the functions
:func:`rfft2` and :func:`irfft2` for 2-D real transforms;
:func:`rfftn` and :func:`irfftn` for N-D real transforms.

The example below demonstrates a 2-D IFFT and plots the resulting
(2-D) time-domain signals.

.. plot::
    :alt: "This code generates six heatmaps arranged in a 2x3 grid. The top row shows mostly blank canvases with the exception of two tiny red peaks on each image. The bottom row shows the real-part of the inverse FFT of each image above it. The first column has two dots arranged horizontally in the top image and in the bottom image a smooth grayscale plot of 5 black vertical stripes representing the 2-D time domain signal. The second column has two dots arranged vertically in the top image and in the bottom image a smooth grayscale plot of 5 horizontal black stripes representing the 2-D time domain signal. In the last column the top image has two dots diagonally located; the corresponding image below has perhaps 20 black stripes at a 60 degree angle."

    >>> from scipy.fft import ifftn
    >>> import matplotlib.pyplot as plt
    >>> import matplotlib.cm as cm
    >>> import numpy as np
    >>> N = 30
    >>> f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    >>> xf = np.zeros((N,N))
    >>> xf[0, 5] = 1
    >>> xf[0, N-5] = 1
    >>> Z = ifftn(xf)
    >>> ax1.imshow(xf, cmap=cm.Reds)
    >>> ax4.imshow(np.real(Z), cmap=cm.gray)
    >>> xf = np.zeros((N, N))
    >>> xf[5, 0] = 1
    >>> xf[N-5, 0] = 1
    >>> Z = ifftn(xf)
    >>> ax2.imshow(xf, cmap=cm.Reds)
    >>> ax5.imshow(np.real(Z), cmap=cm.gray)
    >>> xf = np.zeros((N, N))
    >>> xf[5, 10] = 1
    >>> xf[N-5, N-10] = 1
    >>> Z = ifftn(xf)
    >>> ax3.imshow(xf, cmap=cm.Reds)
    >>> ax6.imshow(np.real(Z), cmap=cm.gray)
    >>> plt.show()


Discrete Cosine Transforms
--------------------------

SciPy provides a DCT with the function :func:`dct` and a corresponding IDCT
with the function :func:`idct`. There are 8 types of the DCT [WPC]_, [Mak]_;
however, only the first 4 types are implemented in scipy. "The" DCT generally
refers to DCT type 2, and "the" Inverse DCT generally refers to DCT type 3. In
addition, the DCT coefficients can be normalized differently (for most types,
scipy provides ``None`` and ``ortho``). Two parameters of the dct/idct
function calls allow setting the DCT type and coefficient normalization.

For a single dimension array x, dct(x, norm='ortho') is equal to
MATLAB dct(x).


Type I DCT
__________

SciPy uses the following definition of the unnormalized DCT-I
(``norm=None``):

.. math::

    y[k] = x_0 + (-1)^k x_{N-1} + 2\sum_{n=1}^{N-2} x[n]
    \cos\left(\frac{\pi nk}{N-1}\right),
    \qquad 0 \le k < N.

Note that the DCT-I is only supported for input size > 1.

Type II DCT
___________

SciPy uses the following definition of the unnormalized DCT-II
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \cos \left({\pi(2n+1)k \over 2N} \right)
    \qquad 0 \le k < N.

In case of the normalized DCT (``norm='ortho'``), the DCT coefficients
:math:`y[k]` are multiplied by a scaling factor `f`:

.. math::

    f = \begin{cases} \sqrt{1/(4N)}, & \text{if $k = 0$} \\    \sqrt{1/(2N)},
    & \text{otherwise} \end{cases} \, .

In this case, the DCT "base functions" :math:`\phi_k[n] = 2 f \cos
\left({\pi(2n+1)k \over 2N} \right)` become orthonormal:

.. math::

   \sum_{n=0}^{N-1} \phi_k[n] \phi_l[n] = \delta_{lk}.


Type III DCT
____________

SciPy uses the following definition of the unnormalized DCT-III
(``norm=None``):

.. math::

    y[k] = x_0 + 2 \sum_{n=1}^{N-1} x[n] \cos\left({\pi n(2k+1) \over 2N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = {x_0\over\sqrt{N}} + {2\over\sqrt{N}} \sum_{n=1}^{N-1} x[n]
    \cos\left({\pi n(2k+1) \over 2N}\right) \qquad 0 \le k < N.


Type IV DCT
___________

SciPy uses the following definition of the unnormalized DCT-IV
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \cos\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = \sqrt{2\over N}\sum_{n=0}^{N-1} x[n] \cos\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N


DCT and IDCT
____________


The (unnormalized) DCT-III is the inverse of the (unnormalized) DCT-II, up to a
factor of ``2N``. The orthonormalized DCT-III is exactly the inverse of the
orthonormalized DCT- II. The function :func:`idct` performs the mappings between
the DCT and IDCT types, as well as the correct normalization.

The following example shows the relation between DCT and IDCT for different
types and normalizations.

>>> from scipy.fft import dct, idct
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])

The DCT-II and DCT-III are each other's inverses, so for an orthonormal transform
we return back to the original signal.

>>> dct(dct(x, type=2, norm='ortho'), type=3, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Doing the same under default normalization, however, we pick up an extra scaling
factor of :math:`2N=10` since the forward transform is unnormalized.

>>> dct(dct(x, type=2), type=3)
array([ 10.,  20.,  10., -10.,  15.])

For this reason, we should use the function `idct` using the same type for both,
giving a correctly normalized result.

>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=2), type=2)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Analogous results can be seen for the DCT-I, which is its own inverse up to a
factor of :math:`2(N-1).`

>>> dct(dct(x, type=1, norm='ortho'), type=1, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>> # Unnormalized round-trip via DCT-I: scaling factor 2*(N-1) = 8
>>> dct(dct(x, type=1), type=1)
array([ 8. ,  16.,  8. , -8. ,  12.])
>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=1), type=1)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

And for the DCT-IV, which is also its own inverse up to a factor of :math:`2N.`

>>> dct(dct(x, type=4, norm='ortho'), type=4, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>> # Unnormalized round-trip via DCT-IV: scaling factor 2*N = 10
>>> dct(dct(x, type=4), type=4)
array([ 10.,  20.,  10., -10.,  15.])
>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=4), type=4)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Example
_______

The DCT exhibits the "energy compaction property", meaning that for many
signals only the first few DCT coefficients have significant magnitude.
Zeroing out the other coefficients leads to a small reconstruction error, a
fact which is exploited in lossy signal compression (e.g. JPEG compression).

The example below shows a signal x and two reconstructions (:math:`x_{20}` and
:math:`x_{15}`) from the signal's DCT coefficients. The signal :math:`x_{20}`
is reconstructed from the first 20 DCT coefficients, :math:`x_{15}` is
reconstructed from the first 15 DCT coefficients. It can be seen that the
relative error of using 20 coefficients is still very small (~0.1%), but
provides a five-fold compression rate.


.. plot::
    :alt: "This code generates an X-Y plot showing amplitude on the Y axis and time on the X axis. The first blue trace is the original signal and starts at amplitude 1 and oscillates down to 0 amplitude over the duration of the plot resembling a frequency chirp. The second red trace is the x_20 reconstruction using the DCT and closely follows the original signal in the high amplitude region but it is unclear to the right side of the plot. The third green trace is the x_15 reconstruction using the DCT and is less precise than the x_20 reconstruction but still similar to x."

    >>> from scipy.fft import dct, idct
    >>> import matplotlib.pyplot as plt
    >>> N = 100
    >>> t = np.linspace(0,20,N, endpoint=False)
    >>> x = np.exp(-t/3)*np.cos(2*t)
    >>> y = dct(x, norm='ortho')
    >>> window = np.zeros(N)
    >>> window[:20] = 1
    >>> yr = idct(y*window, norm='ortho')
    >>> sum(abs(x-yr)**2) / sum(abs(x)**2)
    0.0009872817275276098
    >>> plt.plot(t, x, '-bx')
    >>> plt.plot(t, yr, 'ro')
    >>> window = np.zeros(N)
    >>> window[:15] = 1
    >>> yr = idct(y*window, norm='ortho')
    >>> sum(abs(x-yr)**2) / sum(abs(x)**2)
    0.06196643004256714
    >>> plt.plot(t, yr, 'g+')
    >>> plt.legend(['x', '$x_{20}$', '$x_{15}$'])
    >>> plt.grid()
    >>> plt.show()

Discrete Sine Transforms
------------------------

SciPy provides a DST [Mak]_ with the function :func:`dst` and a corresponding IDST
with the function :func:`idst`.

There are, theoretically, 8 types of the DST for different combinations of
even/odd boundary conditions and boundary offsets [WPS]_, only the first 4
types are implemented in scipy.

Type I DST
__________

DST-I assumes the input is odd around n=-1 and n=N. SciPy uses the following
definition of the unnormalized DST-I (``norm=None``):

.. math::

    y[k] = 2\sum_{n=0}^{N-1} x[n]  \sin\left( \pi {(n+1) (k+1)}\over{N+1}
    \right), \qquad 0 \le k < N.

Note also that the DST-I is only supported for input size > 1. The
(unnormalized) DST-I is its own inverse, up to a factor of ``2(N+1)``.

Type II DST
___________

DST-II assumes the input is odd around n=-1/2 and even around n=N. SciPy uses
the following definition of the unnormalized DST-II (``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n]  \sin\left( {\pi (n+1/2)(k+1)} \over N
    \right), \qquad 0 \le k < N.

Type III DST
____________

DST-III assumes the input is odd around n=-1 and even around n=N-1. SciPy uses
the following definition of the unnormalized DST-III (``norm=None``):

.. math::

    y[k] = (-1)^k x[N-1] + 2 \sum_{n=0}^{N-2} x[n] \sin \left( {\pi
    (n+1)(k+1/2)} \over N \right), \qquad 0 \le k < N.

Type IV DST
___________

SciPy uses the following definition of the unnormalized DST-IV
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \sin\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = \sqrt{2\over N}\sum_{n=0}^{N-1} x[n] \sin\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,


DST and IDST
____________


The following example shows the relation between DST and IDST for
different types and normalizations.

>>> from scipy.fft import dst, idst
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])

The DST-II and DST-III are each other's inverses, so for an orthonormal transform
we return back to the original signal.

>>> dst(dst(x, type=2, norm='ortho'), type=3, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Doing the same under default normalization, however, we pick up an extra scaling
factor of :math:`2N=10` since the forward transform is unnormalized.

>>> dst(dst(x, type=2), type=3)
array([ 10.,  20.,  10., -10.,  15.])

For this reason, we should use the function `idst` using the same type for both,
giving a correctly normalized result.

>>> idst(dst(x, type=2), type=2)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Analogous results can be seen for the DST-I, which is its own inverse up to a
factor of :math:`2(N-1).`

>>> dst(dst(x, type=1, norm='ortho'), type=1, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>>  # scaling factor 2*(N+1) = 12
>>> dst(dst(x, type=1), type=1)
array([ 12.,  24.,  12., -12.,  18.])
>>>  # no scaling factor
>>> idst(dst(x, type=1), type=1)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

And for the DST-IV, which is also its own inverse up to a factor of :math:`2N.`

>>> dst(dst(x, type=4, norm='ortho'), type=4, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>>  # scaling factor 2*N = 10
>>> dst(dst(x, type=4), type=4)
array([ 10.,  20.,  10., -10.,  15.])
>>>  # no scaling factor
>>> idst(dst(x, type=4), type=4)
array([ 1. ,  2. ,  1. , -1. ,  1.5])


Fast Hankel Transform
---------------------

SciPy provides the functions ``fht`` and ``ifht`` to perform the Fast
Hankel Transform (FHT) and its inverse (IFHT) on logarithmically-spaced input
arrays.

The FHT is the discretised version of the continuous Hankel transform defined
by [Ham00]_

.. math::

    A(k) = \int_{0}^{\infty} \! a(r) \, J_{\mu}(kr) \, k \, dr \;,

with :math:`J_{\mu}` the Bessel function of order :math:`\mu`. Under a change
of variables :math:`r \to \log r`, :math:`k \to \log k`, this becomes

.. math::

    A(e^{\log k})
    = \int_{0}^{\infty} \! a(e^{\log r}) \, J_{\mu}(e^{\log k + \log r})
                                        \, e^{\log k + \log r} \, d{\log r}

which is a convolution in logarithmic space. The FHT algorithm uses the FFT
to perform this convolution on discrete input data.

Care must be taken to minimise numerical ringing due to the circular nature
of FFT convolution. To ensure that the low-ringing condition [Ham00]_ holds,
the output array can be slightly shifted by an offset computed using the
``fhtoffset`` function.


Appendix
--------

The applications of Fourier transforms are too many to be listed here.
The aim of this appendix is to give some insight on how to interpret
a spectrum, rather than to introduce domain specific applications.

Spectral leakage
________________

Let :math:`x[n]` be a sequence of length :math:`N,` which represents samples
from a continuous-time signal :math:`x(t),` taken at time :math:`t=n\,T_s.`
Its DFT :math:`y[k]` represents the complex spectum of these samples at
frequency :math:`f=k\,F_p,` where :math:`F_p = 1/(N\,T_s)` is the bin width.
The sequence :math:`y[k]` can be extended to a continuous-frequency function
:math:`y(f),` by substituting :math:`k=f/F_p` in the definition of DFT.
This function is periodic with period equal to :math:`F_s=1/T_s.`
We can consider that :math:`y(f)` reveals the complex spectrum at
fractional values of :math:`k=f/F_p.` Notice that using the inverse DFT,
the function :math:`y(f)` can be expressed in terms of the sequence
:math:`y[k],` therefore both contain the same information about the
sequence :math:`x[n].`

An easy way to calculate the complex spectrum :math:`y(f),` is to append
:math:`(M-1)\,N` zeros at the end of the sequence :math:`x[n].`
The extended sequence has length :math:`M\,N` and its FFT gives the
values of :math:`y(f)` with resolution :math:`F_p/M,`
as shown in the following example.

.. plot::
    :alt: "This code calculates the spectrum of a signal at high resolution, using zero-extension. and plots its amplitude in dB, i.e., in 20*log10. The spectrum contains three peaks, at zero frequency, at 50 Hz and at 125 Hz. Two inset axes show more detail around the peaks at 50 Hz and 125 Hz."

    >>> from scipy.fft import rfft, rfftfreq
    >>> import numpy as np
    >>> M = 10         # number of segments
    >>> N = 400        # number of samples per segment
    >>> Fs = 800.0     # sampling frequency in Hz
    >>> Ts = 1.0 / Fs  # sampling time (step) in s
    >>> Tp = N / Fs    # sampling period in s (one segment)
    >>> Fp = Fs / N    # main bin width in Hz
    >>> Fq = Fp / M    # sub-bin width in Hz
    >>> Fx = np.array([0.0, 50.0, 125.0])            # signal frequencies in Hz
    >>> Ax = np.array([0.2, 0.8, 0.6])               # signal amplitudes
    >>> t = np.linspace(0.0, Tp, N, endpoint=False)  # time in s (one segment)
    >>> f = rfftfreq(M*N, Ts)                        # frequency in Hz
    >>> x = np.zeros(M*N)
    >>> x[:N] = np.sum(Ax[:,np.newaxis]*np.cos(2.0*np.pi*np.outer(Fx,t)), axis=0)
    >>> yr = rfft(x)
    >>> ya = 20.0*np.log10(2.0/N * np.abs(yr))
    >>>
    >>> # plot ya vs. f in the range [0, Fs/2]
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(7,4))
    >>> fig.tight_layout()
    >>> ax.plot(f, ya)
    >>>
    >>> # decorations
    >>> for i in range(len(Fx)):
    ...     a = 20*np.log10(Ax[i] if i > 0 else 2*Ax[i])
    ...     ax.plot(Fx[i]+np.array([[-6.0, 6.0], [-3.0, 3.0]]), a*np.r_[1.0, 1.0], 'g')
    >>> ax.set_ylim(-60.0, 0.0)
    >>> ax.set_title(f'M={M}, N={N}, Fp/M={Fq}, Fp={Fp} Hz, Fs={Fs} Hz')
    >>> ax.set_xlabel('f [Hz]')
    >>> ax.set_ylabel('2/N |y| [dB]')
    >>> ax.grid()
    >>> for i in (1, 2):
    ...     axi = fig.add_axes([0.35+0.22*i, 0.45, 0.13, 0.4])
    ...     k = np.r_[int(Fx[i]/Fq+0.3)-4*M:int(Fx[i]/Fq+0.8)+4*M+1]
    ...     axi.plot(f[k], ya[k])
    ...     axi.set_xlim(f[k[0]], f[k[-1]])
    ...     axi.set_ylim(-60.0, 0.0)
    ...     axi.grid()
    >>> plt.show()

The amplitude of the spectrum is shown in dB, i.e., in 20*log10.
As it is revealed by this example, each sinusoidal component in the signal,
including the bias at zero frequency, appears with a common pattern around
its frequency. This pattern is called "spectral leakage".

For a complex exponential signal of the form

.. math::

    x[n] = e^{j 2 \pi F_x n T_s} \, ,

where :math:`F_x` is the signal frequency, it can be shown that

.. math::

    \big|\frac{1}{N}\,y(f)\big| = |{\rm sinc}_N(\frac{f-F_x}{F_p})| \, ,

where

.. math::

    {\rm sinc}_N(u) = \frac{\sin(\pi u)}{N\,\sin(\pi u/N)}

is the "periodic sampling function" with period :math:`N.`

This pattern has a main peak at the signal frequency :math:`F_x,`
and zeros at all integer multiples of :math:`F_p` relative to :math:`F_x,`
excluding integer multiples of :math:`F_s`
(since the spectrum has period :math:`F_s`).
In the worst case, the signal frequency is in the middle between
two bin frequencies, and the amplitude in these two bins is scaled
by :math:`2/\pi`, i.e., it appears 4 dB below the actual peak.
The side peaks are close to integer-plus-half multiples of :math:`F_p`
relative to :math:`F_x.`
When the spectrum is plotted in dB, the envelop of the spectral leakage
near the main peak takes the form of a negative logarithmic function,
which falls by 6 dB (or 1 bit in log2) per doubling of the distance
from the main peak.
The second peak is 13 dB below the main peak, while the spectral leakage 
at the far end (at distance :math:`F_s/2`) is :math:`1/N.`
This spectral leakage is prohibitively high for most signal processing
applications.

For :math:`F_x=0,` the signal is :math:`x[n]=1,` therefore the
spectral leakage pattern is the spectrum of the rectangular window
(:math:`N` ones followed by zeros).


Windows
_______

The standard way to reduce spectral leakage is to pre-mulitply the signal
sequence :math:`x[n]` with a "window" sequence :math:`w[n]` of the same length,
before the calculation of FFT.
Typically the window sequence :math:`w[n]` is a smooth function of :math:`n,`
with a peak at the center and zero or small values at the ends of the range
of :math:`n,` such that the transition into and out of the sampling period
is much smoother that with the (implied) rectangular window.
There is a big number of useful windows, with different trade-offs between
their performance parameters (see [WPW]_ and [Har78]_).

Two of the windows available in scipy.signal.windows are shown below as an example:
the Blackman window and the Dolph–Chebyshev window with 100 dB attenuation.

.. plot::
    :alt: "This code generates a plot of (1) the  Blackman window, and (2) the Dolph–Chebyshev window with 100 dB attenuation. Both window functions have a peak at the center of the horizontal axis. Their maximum value is normalized to 1."

    >>> from scipy.signal import windows
    >>> import numpy as np
    >>> N = 400  # window length
    >>> n = np.r_[:N]
    >>> # plot w[n] vs. n
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(7,4))
    >>> fig.tight_layout()
    >>> ax.plot(n, windows.blackman(N), 'g', label='blackman')
    >>> ax.plot(n, windows.chebwin(N, at=100), 'r', label='cheb(100)')
    >>> ax.legend()
    >>>
    >>> # decorations
    >>> ax.set_title(f'N={N}')
    >>> ax.set_xlabel('n')
    >>> ax.set_ylabel('w[n]')
    >>> ax.grid()
    >>> plt.show()

Their spectrum, which shows also their spectral leakage pattern,
is calculated and plotted with the following code:

.. plot::
    :alt: "This code generates calculates and plots the spectrum of three windows: (1) the rectangular window, (2) the  Blackman window, and (3) the Dolph–Chebyshev window with 100 dB attenuation. The main axes shows the spectrum at integer-plus-half bin width, which is approximately equal to the envelop of the spectral leakage. An inset axes shows a detail of the spectrum for the first 6 bins."

    >>> from scipy.fft import rfft, rfftfreq
    >>> from scipy.signal import windows
    >>> import numpy as np
    >>> M = 10   # number of segments
    >>> N = 400  # window length
    >>> f = rfftfreq(M*N, 1.0/N)  # normalized frequency
    >>> x = np.zeros((3, M*N))    # space for 3 windows
    >>> x[0,:N] = np.ones(N)
    >>> x[1,:N] = windows.blackman(N)
    >>> x[2,:N] = windows.chebwin(N, at=100)
    >>> yr = rfft(x, axis=1)
    >>> eps = np.spacing(1.0)
    >>> ya = 20.0*np.log10(eps + np.abs(yr) / np.sum(x, axis=1, keepdims=True))
    >>> # envelop of ya
    >>> fm = f[:-1].reshape(-1,M)
    >>> yam = ya[:,:-1].reshape(3,-1,M)
    >>> ke = np.argmax(yam, axis=-1)
    >>> fe = np.take_along_axis(fm[np.newaxis], ke[:,:,np.newaxis], axis=-1).squeeze(axis=-1)
    >>> yae = np.take_along_axis(yam, ke[:,:,np.newaxis], axis=-1).squeeze(axis=-1)
    >>>
    >>> # plot normalized spectrum vs. normalized frequency
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(7,4))
    >>> fig.tight_layout()
    >>> ax.set_prop_cycle(color=['b', 'g', 'r'])
    >>> ax.plot(fe.T, yae.T)
    >>> ax.legend(['rectangle', 'blackman', 'cheb(100)'])
    >>>
    >>> # decorations
    >>> ax.set_xlim(-1, N//2+1)
    >>> ax.set_ylim(-150.0, 10.0)
    >>> ax.set_title(f'N={N}, Fp=1, Fs={N}')
    >>> ax.set_xlabel('normalized frequency (f/Fp)')
    >>> ax.set_ylabel('normalized spectrum [dB]')
    >>> ax.grid()
    >>> axi = fig.add_axes([0.5, 0.2, 0.4, 0.5])
    >>> axi.set_prop_cycle(color=['b', 'g', 'r'])
    >>> ki = np.r_[0:M*6+1]
    >>> axi.plot(f[ki], ya[:,ki].T)
    >>> axi.set_xlim(0.0, 6.0)
    >>> axi.set_ylim(-130.0, 10.0)
    >>> axi.grid()
    >>> plt.show()

In the previous plot, the frequency is normalized such that the
bin frequencies are the integer values in the horizontal axis
(this is achieved by setting the sampling time equal to :math:`1/N`).
The amplitude is normalized such that its value at zero frequency is 0 dB
(this is achieved by dividing either the window function or its spectum
by the sum of the window sequence).
The inset shows the amplitude for the first few bins.
As it is seen, the amplitude has large variation over the continuous frequency,
due to the zero values at or near every bin frequency.
The main plot shows the envelop of the amplitude vs. frequency,
which is also the worst-case pattern for the spectral leakage.

These three windows are applied to a signal as shown below.

.. plot::
    :alt: "This code calculates and plots the spectrum of a signal using three different windows: (1) rectangular, (2) Blackman, and (3) Dolph–Chebyshev window with 100 dB attenuation. The spectrum contains three peaks: at zero frequency, at 50 Hz and at 125 Hz. An inset axes shows the spectral leakage in more detail, around the peak at 125 Hz."

    >>> from scipy.fft import rfft, rfftfreq
    >>> from scipy.signal import windows
    >>> import numpy as np
    >>> N = 400        # number of samples
    >>> Fs = 800.0     # sampling frequency in Hz
    >>> Ts = 1.0 / Fs  # sampling time (step) in s
    >>> Tp = N / Fs    # sampling period (total) in s
    >>> Fp = Fs / N    # bin width (frequency resolution) in Hz
    >>> Fx = np.array([0.0, 50.0, 125.0])            # signal frequencies in Hz
    >>> Ax = np.array([0.2, 0.8, 0.6])               # signal amplitudes
    >>> t = np.linspace(0.0, Tp, N, endpoint=False)  # time in s
    >>> f = rfftfreq(N, Ts)                          # frequency in Hz
    >>> w = np.empty((3, N))                         # space for 3 windows
    >>> w[0,:N] = np.ones(N)
    >>> w[1,:N] = windows.blackman(N)
    >>> w[2,:N] = windows.chebwin(N, at=100)
    >>> x = w * np.sum(Ax[:,np.newaxis]*np.cos(2.0*np.pi*np.outer(Fx,t)), axis=0)
    >>> yr = rfft(x, axis=1)
    >>> eps = np.spacing(1.0)
    >>> ya = 20.0*np.log10(eps + 2.0 * np.abs(yr) / np.sum(w, axis=1, keepdims=True))
    >>>
    >>> # plot normalized spectrum vs. f for different windows
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(figsize=(7,4))
    >>> fig.tight_layout()
    >>> ax.set_prop_cycle(color=['b', 'g', 'r'])
    >>> ax.plot(f, ya.T)
    >>> ax.legend(['rectangle', 'blackman', 'cheb(100)'])
    >>>
    >>> # decorations
    >>> for i in range(len(Fx)):
    ...     a = 20.0*np.log10(Ax[i] if i > 0 else 2*Ax[i])
    ...     ax.plot(Fx[i]+np.array([[-6.0, 6.0], [-3.0, 3.0]]), a*np.r_[1.0, 1.0], 'g')
    >>> ax.set_ylim(-120.0, 0.0)
    >>> ax.set_title(f'N={N}, Fp={Fp} Hz, Fs={Fs} Hz')
    >>> ax.set_xlabel('f [Hz]')
    >>> ax.set_ylabel('amplitude [dB]')
    >>> ax.grid()
    >>> axi = fig.add_axes([0.55, 0.22, 0.3, 0.5])
    >>> axi.set_prop_cycle(color=['b', 'g', 'r'])
    >>> ki = np.r_[int(Fx[2]/Fp+0.3)-10:int(Fx[2]/Fp+0.8)+11]
    >>> axi.plot(f[ki], ya[:,ki].T)
    >>> axi.set_xlim(f[ki[0]], f[ki[-1]])
    >>> axi.set_ylim(-120.0, 0.0)
    >>> axi.grid()
    >>> plt.show()

As it is seen in the last example, if the frequency of a sinusoidal
component is an exact multiple of the bin width, then this component
appears as an impulse when a rectangular window (no window) is used,
which is the best possible pattern, while with other windows
it appears to be broader than an impulse.
However, except for the bias which is always at zero frequency,
this condition is seldom true (usually also the assumption that the
energy of a signal component is concentated at one frequency is not true).

The main requirements in the design of a window function are the following:

* It must be time-limited, i.e., its duration is equal to a parameter :math:`N.`

* The gain within the first bin must be flat, otherwise the reading of the
  amplitude at a peak does not give a good estimate of the actual (hidden) peak.

* The passband bandwidth must be small, otherwise the width of a sinusoidal
  or narrow-bandwidth component is too big.

* The attenuation outside of the passband must increase rapidly with frequency,
  otherwise the spectral leakage is too big.

In addition, a window function must be easy to calculate for any length :math:`N.`
As it is seen in previous examples, in order to fulfill these requirements,
the window functions are smooth. However, a function which is smooth and
has a bell-like shape, is not necessarily a good window function.
For example, the product of two (same or different) window functions is
also a smooth function, but its passband bandwidth is increased without
much (or any) improvement in stopband attenuation.
This means that it is not a good idea to apply a second window in order
to reduce the spectral leakage.


References
----------

.. [CT65] Cooley, James W., and John W. Tukey, 1965, "An algorithm for the
        machine calculation of complex Fourier series," *Math. Comput.*
        19: 297-301.

.. [NR07] Press, W., Teukolsky, S., Vetterline, W.T., and Flannery, B.P.,
        2007, *Numerical Recipes: The Art of Scientific Computing*, ch.
        12-13.  Cambridge Univ. Press, Cambridge, UK.

.. [Mak] J. Makhoul, 1980, "A fast cosine transform in one and two dimensions",
       `IEEE Trans. Acoustics, Speech, and Signal Processing`,
       vol. 28(1), pp. 27-34, :doi:`10.1109/TASSP.1980.1163351`

.. [Ham00] A. J. S. Hamilton, 2000, "Uncorrelated modes of the non-linear power
       spectrum", *MNRAS*, 312, 257. :doi:`10.1046/j.1365-8711.2000.03071.x`

.. [Har78] F. J. Harris, 1978, "On the use of windows for harmonic analysis
       with the discrete Fourier transform", `Proc. IEEE`, vol. 66(1), pp. 51-83.

.. [WPW] https://en.wikipedia.org/wiki/Window_function

.. [WPC] https://en.wikipedia.org/wiki/Discrete_cosine_transform

.. [WPS] https://en.wikipedia.org/wiki/Discrete_sine_transform
