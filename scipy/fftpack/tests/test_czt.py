# This program is public domain
# Authors: Paul Kienzle, Nadav Horesh
'''
  A unit test module for czt.py
'''

from numpy.testing import run_module_suite
from czt import czt, zoomfft, scaledfft
import numpy as np

fft = np.fft.fft
fftshift = np.fft.fftshift
norm = np.linalg.norm


def check_zoomfft(x):

    # Normal fft and zero-padded fft equivalent to 10x oversampling
    over = 10
    w = np.linspace(0, 2-2./len(x), len(x))
    y = fft(x)
    wover = np.linspace(0, 2-2./(over*len(x)), over*len(x))
    yover = fft(x, over*len(x))

    # Check that zoomfft is the equivalent of fft
    y1 = zoomfft(x, [0, 2-2./len(y)])

    # Check that zoomfft with oversampling is equivalent to zero padding
    y2 = zoomfft(x, [0, 2-2./len(yover)], m=len(yover))

    # Check that zoomfft works on a subrange
    f1, f2 = w[3], w[6]
    y3 = zoomfft(x, [f1, f2], m=3*over+1)
    w3 = np.linspace(f1, f2, len(y3))
    idx3 = slice(3*over, 6*over+1)

    err = norm(y-y1)/norm(y)
    # print "direct err %g"%err
    assert err < 1e-10, "error for direct transform is %g" % (err,)
    err = norm(yover-y2)/norm(yover)
    # print "over err %g"%err
    assert err < 1e-10, "error for oversampling is %g" % (err,)
    err = norm(yover[idx3]-y3)/norm(yover[idx3])
    # print "range err %g"%err
    assert err < 1e-10, "error for subrange is %g" % (err,)


def check_scaledfft(x):
    n = len(x)
    assert norm(fft(x)-scaledfft(x)) < 2E-10
    n = n//4 * 4
    x = x[:n]
    st = n//4
    m = n//2
    sliced_fft = fftshift(fft(x))[st:st+m]
    scl_fft = fftshift(scaledfft(x, scale=0.5, m=m))
    assert norm(sliced_fft-scl_fft) < 2E-10


def test_1D():
    'Test of 1D version of the transforms'
    lengths = np.random.randint(8, 200, 20)
    for length in lengths:
        x = np.random.random(length)
        yield check_zoomfft, x
        yield check_scaledfft, x

if __name__ == '__main__':
    np.random.seed()
    run_module_suite()
