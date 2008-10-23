import numpy as np

"""This module defines the function frame() which creates
a framed copy of an input array with the boundary pixels
defined according to the IRAF boundary modes: 'nearest',
'reflect', 'wrap', and 'constant.'
"""

def frame_nearest(a, shape, cval=None):

    """frame_nearest creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    copied from the nearest edge pixel in 'a'.

    >>> a = np.arange(16)
    >>> a.shape=(4,4)
    >>> frame_nearest(a, (8,8))
    array([[ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 0,  0,  0,  1,  2,  3,  3,  3],
           [ 4,  4,  4,  5,  6,  7,  7,  7],
           [ 8,  8,  8,  9, 10, 11, 11, 11],
           [12, 12, 12, 13, 14, 15, 15, 15],
           [12, 12, 12, 13, 14, 15, 15, 15],
           [12, 12, 12, 13, 14, 15, 15, 15]])

    """

    b = np.zeros(shape, dtype=a.dtype)
    delta = (np.array(b.shape) - np.array(a.shape))
    dy = delta[0] // 2
    dx = delta[1] // 2
    my = a.shape[0] + dy
    mx = a.shape[1] + dx

    b[dy:my, dx:mx] = a                  # center
    b[:dy,dx:mx]  = a[0:1,:]               # top
    b[my:,dx:mx]  = a[-1:,:]              # bottom
    b[dy:my, :dx] = a[:, 0:1]              # left
    b[dy:my, mx:] = a[:, -1:]             # right
    b[:dy, :dx]   = a[0,0]               # topleft
    b[:dy, mx:]   = a[0,-1]              # topright
    b[my:, :dx]   = a[-1, 0]             # bottomleft
    b[my:, mx:]   = a[-1, -1]            # bottomright

    return b

def frame_reflect(a, shape, cval=None):

    """frame_reflect creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    reflected from the nearest edge pixels in 'a'.

    >>> a = np.arange(16)
    >>> a.shape = (4,4)
    >>> frame_reflect(a, (8,8))
    array([[ 5,  4,  4,  5,  6,  7,  7,  6],
           [ 1,  0,  0,  1,  2,  3,  3,  2],
           [ 1,  0,  0,  1,  2,  3,  3,  2],
           [ 5,  4,  4,  5,  6,  7,  7,  6],
           [ 9,  8,  8,  9, 10, 11, 11, 10],
           [13, 12, 12, 13, 14, 15, 15, 14],
           [13, 12, 12, 13, 14, 15, 15, 14],
           [ 9,  8,  8,  9, 10, 11, 11, 10]])
    """

    b = np.zeros(shape, dtype=a.dtype)
    delta = (np.array(b.shape) - np.array(a.shape))
    dy = delta[0] // 2
    dx = delta[1] // 2
    my = a.shape[0] + dy
    mx = a.shape[1] + dx
    sy = delta[0] - dy
    sx = delta[1] - dx

    b[dy:my, dx:mx] = a                            # center
    b[:dy,dx:mx]  = a[:dy,:][::-1,:]               # top
    b[my:,dx:mx]  = a[-sy:,:][::-1,:]              # bottom
    b[dy:my,:dx]  = a[:,:dx][:,::-1]               # left
    b[dy:my,mx:]  = a[:,-sx:][:,::-1]              # right
    b[:dy,:dx]    = a[:dy,:dx][::-1,::-1]          # topleft
    b[:dy,mx:]    = a[:dy,-sx:][::-1,::-1]         # topright
    b[my:,:dx]    = a[-sy:,:dx][::-1,::-1]         # bottomleft
    b[my:,mx:]    = a[-sy:,-sx:][::-1,::-1]        # bottomright
    return b

def frame_wrap(a, shape, cval=None):
    """frame_wrap creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    wrapped around to the opposite edge pixels in 'a'.

    >>> a = np.arange(16)
    >>> a.shape=(4,4)
    >>> frame_wrap(a, (8,8))
    array([[10, 11,  8,  9, 10, 11,  8,  9],
           [14, 15, 12, 13, 14, 15, 12, 13],
           [ 2,  3,  0,  1,  2,  3,  0,  1],
           [ 6,  7,  4,  5,  6,  7,  4,  5],
           [10, 11,  8,  9, 10, 11,  8,  9],
           [14, 15, 12, 13, 14, 15, 12, 13],
           [ 2,  3,  0,  1,  2,  3,  0,  1],
           [ 6,  7,  4,  5,  6,  7,  4,  5]])

    """

    b = np.zeros(shape, dtype=a.dtype)
    delta = (np.array(b.shape) - np.array(a.shape))
    dy = delta[0] // 2
    dx = delta[1] // 2
    my = a.shape[0] + dy
    mx = a.shape[1] + dx
    sy = delta[0] - dy
    sx = delta[1] - dx

    b[dy:my, dx:mx] = a                  # center
    b[:dy,dx:mx]  = a[-dy:,:]            # top
    b[my:,dx:mx]  = a[:sy,:]             # bottom
    b[dy:my,:dx]  = a[:,-dx:]            # left
    b[dy:my,mx:]  = a[:, :sx]            # right
    b[:dy,:dx]    = a[-dy:,-dx:]         # topleft
    b[:dy,mx:]    = a[-dy:,:sx ]         # topright
    b[my:,:dx]    = a[:sy, -dx:]         # bottomleft
    b[my:,mx:]    = a[:sy, :sx]          # bottomright
    return b

def frame_constant(a, shape, cval=0):
    """frame_nearest creates an oversized copy of 'a' with new 'shape'
    and the contents of 'a' in the center.  The boundary pixels are
    copied from the nearest edge pixel in 'a'.

    >>> a = np.arange(16)
    >>> a.shape=(4,4)
    >>> frame_constant(a, (8,8), cval=42)
    array([[42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42,  0,  1,  2,  3, 42, 42],
           [42, 42,  4,  5,  6,  7, 42, 42],
           [42, 42,  8,  9, 10, 11, 42, 42],
           [42, 42, 12, 13, 14, 15, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42],
           [42, 42, 42, 42, 42, 42, 42, 42]])

    """

    b = np.zeros(shape, dtype=a.dtype)
    delta = (np.array(b.shape) - np.array(a.shape))
    dy = delta[0] // 2
    dx = delta[1] // 2
    my = a.shape[0] + dy
    mx = a.shape[1] + dx

    b[dy:my, dx:mx] = a              # center
    b[:dy,dx:mx]  = cval             # top
    b[my:,dx:mx]  = cval             # bottom
    b[dy:my, :dx] = cval             # left
    b[dy:my, mx:] = cval             # right
    b[:dy, :dx]   = cval             # topleft
    b[:dy, mx:]   = cval             # topright
    b[my:, :dx]   = cval             # bottomleft
    b[my:, mx:]   = cval             # bottomright
    return b

_frame_dispatch = { "nearest": frame_nearest,
                    "reflect": frame_reflect,
                    "wrap": frame_wrap,
                    "constant" : frame_constant }

def frame(a, shape, mode="nearest", cval=0.0):

    """frame creates an oversized copy of 'a' with new 'shape', with
    extra pixels being supplied according to IRAF boundary mode,
    'mode'.  """

    try:
        f = _frame_dispatch[mode]
    except KeyError:
        raise ValueError('invalid IRAF boundary mode: "%s"' % mode)

    return f(a, shape, cval)

def unframe(a, shape):

    """unframe extracts the center slice of framed array 'a' which had
    'shape' prior to framing."""

    delta = np.array(a.shape) - np.array(shape)
    dy = delta[0]//2
    dx = delta[1]//2
    my = shape[0] + dy
    mx = shape[1] + dx
    return a[dy:my, dx:mx]

def test():
    import doctest, iraf_frame
    return doctest.testmod(iraf_frame)
