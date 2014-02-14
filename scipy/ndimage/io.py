from __future__ import division, print_function, absolute_import

__all__ = ['imread']

from numpy import array


def imread(fname, flatten=False, mode=None):
    """
    Read an image from a file as an array.

    Parameters
    ----------
    fname : str
        Image file name, e.g. ``test.jpg``, or a file object.
    flatten : bool, optional
        If true, convert the output to grey-scale. Default is False.
    mode : str, optional
        mode to convert image to, e.g. ``RGB``.


    Returns
    -------
    img_array : ndarray
        The different colour bands/channels are stored in the
        third dimension, such that a grey-image is MxN, an
        RGB-image MxNx3 and an RGBA-image MxNx4.

    Raises
    ------
    ImportError
        If the Python Imaging Library (PIL) can not be imported.

    """
    try:
        from PIL import Image
    except ImportError:
        raise ImportError("Could not import the Python Imaging Library (PIL)"
                          " required to load image files.  Please refer to"
                          " http://pypi.python.org/pypi/PIL/ for installation"
                          " instructions.")

    im = Image.open(fname)
    if mode:
        im = im.convert(mode)
    if flatten:
        im = im.convert('F')
    result = array(im)
    return result
