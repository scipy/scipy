from numpy import array

def imread(fname, flatten=False):
    """Load an image from file.

    Parameters
    ----------
    im : PIL image
        Input image.
    flatten : bool
        If true, convert the output to grey-scale.

    Returns
    -------
    img_array : ndarray
        The different colour bands/channels are stored in the
        third dimension, such that a grey-image is MxN, an
        RGB-image MxNx3 and an RGBA-image MxNx4.

    """
    from PIL import Image
    im = Image.open(fname)
    if flatten:
        im = im.convert('F')
    return array(im)

