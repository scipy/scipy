# Functions which need the PIL

import numpy
import tempfile

from numpy import amin, amax, ravel, asarray, cast, arange, \
     ones, newaxis, transpose, mgrid, iscomplexobj, sum, zeros, uint8, \
     issubdtype, array

import Image
import ImageFilter

__all__ = ['fromimage','toimage','imsave','imread','bytescale',
           'imrotate','imresize','imshow','imfilter','radon']

# Returns a byte-scaled image
def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    """
    Parameters
    ----------
    im : PIL image
         Input image.
    flatten : bool
              If true, convert the output to grey-scale

    Returns
    -------
    img_array : ndarray
                The different colour bands/channels are stored in the
                third dimension, such that a grey-image is MxN, an
                RGB-image MxNx3 and an RGBA-image MxNx4.

    """
    if data.dtype == uint8:
        return data
    high = high - low
    if cmin is None: cmin = data.min()
    if cmax is None: cmax = data.max()
    scale = high *1.0 / (cmax-cmin or 1)
    bytedata = ((data*1.0-cmin)*scale + 0.4999).astype(uint8)
    return bytedata + cast[uint8](low)

def imread(name,flatten=0):
    """
    Read an image file from a filename.

    Parameters
    ----------
    name : str
        The file name to be read.
    flatten : bool, optional
        If True, flattens the color layers into a single gray-scale layer.

    Returns
    -------
     : nd_array
        The array obtained by reading image.

    Notes
    -----
    The image is flattened by calling convert('F') on
    the resulting image object.

    """

    im = Image.open(name)
    return fromimage(im,flatten=flatten)

def imsave(name, arr):
    """
    Save an array to an image file.

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
    im = toimage(arr)
    im.save(name)
    return

def fromimage(im, flatten=0):
    """
    Return a copy of a PIL image as a numpy array.

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
    if not Image.isImageType(im):
        raise TypeError("Input is not a PIL image.")
    if flatten:
        im = im.convert('F')
    return array(im)

_errstr = "Mode is unknown or incompatible with input array shape."
def toimage(arr,high=255,low=0,cmin=None,cmax=None,pal=None,
            mode=None,channel_axis=None):
    """Takes a numpy array and returns a PIL image.  The mode of the
    PIL image depends on the array shape, the pal keyword, and the mode
    keyword.

    For 2-D arrays, if pal is a valid (N,3) byte-array giving the RGB values
    (from 0 to 255) then mode='P', otherwise mode='L', unless mode is given
    as 'F' or 'I' in which case a float and/or integer array is made

    For 3-D arrays, the channel_axis argument tells which dimension of the
      array holds the channel data.
    For 3-D arrays if one of the dimensions is 3, the mode is 'RGB'
      by default or 'YCbCr' if selected.
    if the

    The numpy array must be either 2 dimensional or 3 dimensional.
    """
    data = asarray(arr)
    if iscomplexobj(data):
        raise ValueError, "Cannot convert a complex-valued array."
    shape = list(data.shape)
    valid = len(shape)==2 or ((len(shape)==3) and \
                              ((3 in shape) or (4 in shape)))
    assert valid, "Not a suitable array shape for any mode."
    if len(shape) == 2:
        shape = (shape[1],shape[0]) # columns show up first
        if mode == 'F':
            data32 = data.astype(numpy.float32)
            image = Image.fromstring(mode,shape,data32.tostring())
            return image
        if mode in [None, 'L', 'P']:
            bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
            image = Image.fromstring('L',shape,bytedata.tostring())
            if pal is not None:
                image.putpalette(asarray(pal,dtype=uint8).tostring())
                # Becomes a mode='P' automagically.
            elif mode == 'P':  # default gray-scale
                pal = arange(0,256,1,dtype=uint8)[:,newaxis] * \
                      ones((3,),dtype=uint8)[newaxis,:]
                image.putpalette(asarray(pal,dtype=uint8).tostring())
            return image
        if mode == '1':  # high input gives threshold for 1
            bytedata = (data > high)
            image = Image.fromstring('1',shape,bytedata.tostring())
            return image
        if cmin is None:
            cmin = amin(ravel(data))
        if cmax is None:
            cmax = amax(ravel(data))
        data = (data*1.0 - cmin)*(high-low)/(cmax-cmin) + low
        if mode == 'I':
            data32 = data.astype(numpy.uint32)
            image = Image.fromstring(mode,shape,data32.tostring())
        else:
            raise ValueError, _errstr
        return image

    # if here then 3-d array with a 3 or a 4 in the shape length.
    # Check for 3 in datacube shape --- 'RGB' or 'YCbCr'
    if channel_axis is None:
        if (3 in shape):
            ca = numpy.flatnonzero(asarray(shape) == 3)[0]
        else:
            ca = numpy.flatnonzero(asarray(shape) == 4)
            if len(ca):
                ca = ca[0]
            else:
                raise ValueError, "Could not find channel dimension."
    else:
        ca = channel_axis

    numch = shape[ca]
    if numch not in [3,4]:
        raise ValueError, "Channel axis dimension is not valid."

    bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
    if ca == 2:
        strdata = bytedata.tostring()
        shape = (shape[1],shape[0])
    elif ca == 1:
        strdata = transpose(bytedata,(0,2,1)).tostring()
        shape = (shape[2],shape[0])
    elif ca == 0:
        strdata = transpose(bytedata,(1,2,0)).tostring()
        shape = (shape[2],shape[1])
    if mode is None:
        if numch == 3: mode = 'RGB'
        else: mode = 'RGBA'


    if mode not in ['RGB','RGBA','YCbCr','CMYK']:
        raise ValueError, _errstr

    if mode in ['RGB', 'YCbCr']:
        assert numch == 3, "Invalid array shape for mode."
    if mode in ['RGBA', 'CMYK']:
        assert numch == 4, "Invalid array shape for mode."

    # Here we know data and mode is coorect
    image = Image.fromstring(mode, shape, strdata)
    return image

def imrotate(arr,angle,interp='bilinear'):
    """
    Rotate an image counter-clockwise by angle degrees.

    Parameters
    ----------
    arr : nd_array
        Input array of image to be rotated.
    angle : float
        The angle of rotation.
    interp : str, optional
        Interpolation


    Returns
    -------
     : nd_array
        The rotated array of image.

    Notes
    -----

    Interpolation methods can be:
    * 'nearest' :  for nearest neighbor
    * 'bilinear' : for bilinear
    * 'cubic' : cubic
    * 'bicubic' : for bicubic

    """
    arr = asarray(arr)
    func = {'nearest':0,'bilinear':2,'bicubic':3,'cubic':3}
    im = toimage(arr)
    im = im.rotate(angle,resample=func[interp])
    return fromimage(im)

def imresize(arr,newsize,interp='bilinear',mode=None):
    newsize=list(newsize)
    newsize.reverse()
    newsize = tuple(newsize)
    arr = asarray(arr)
    func = {'nearest':0,'bilinear':2,'bicubic':3,'cubic':3}
    im = toimage(arr,mode=mode)
    im = im.resize(newsize,resample=func[interp])
    return fromimage(im)

def imshow(arr):
    """Simple showing of an image through an external viewer.
    """
    im = toimage(arr)
    fnum,fname = tempfile.mkstemp('.png')
    try:
        im.save(fname)
    except:
        raise RuntimeError("Error saving temporary image data.")

    import os
    os.close(fnum)

    cmd = os.environ.get('SCIPY_PIL_IMAGE_VIEWER','see')
    status = os.system("%s %s" % (cmd,fname))

    os.unlink(fname)
    if status != 0:
        raise RuntimeError('Could not execute image viewer.')

def imresize(arr,size):
    """
    Resize an image.

    Parameters
    ----------
    arr : nd_array
        The array of image to be resized.

    size : int, float or tuple
        * int   - Percentage of current size.
        * float - Fraction of current size.
        * tuple - Size of the output image.

    Returns
    -------

     : nd_array
        The resized array of image.

    """
    im = toimage(arr)
    ts = type(size)
    if issubdtype(ts,int):
        size = size / 100.0
    elif issubdtype(type(size),float):
        size = (array(im.size)*size).astype(int)
    else:
        size = (size[1],size[0])
    imnew = im.resize(size)
    return fromimage(imnew)


def imfilter(arr,ftype):
    """
    Simple filtering of an image.

    Parameters
    ----------
    arr : ndarray
        The array of Image in which the filter is to be applied.
    ftype : str
        The filter that has to be applied. Legal values are:
        'blur', 'contour', 'detail', 'edge_enhance', 'edge_enhance_more',
        'emboss', 'find_edges', 'smooth', 'smooth_more', 'sharpen'.

    Returns
    -------
    res : nd_array
        The array with filter applied.

    Raises
    ------
    ValueError
        *Unknown filter type.* . If the filter you are trying
        to apply is unsupported.

    """
    _tdict = {'blur':ImageFilter.BLUR,
              'contour':ImageFilter.CONTOUR,
              'detail':ImageFilter.DETAIL,
              'edge_enhance':ImageFilter.EDGE_ENHANCE,
              'edge_enhance_more':ImageFilter.EDGE_ENHANCE_MORE,
              'emboss':ImageFilter.EMBOSS,
              'find_edges':ImageFilter.FIND_EDGES,
              'smooth':ImageFilter.SMOOTH,
              'smooth_more':ImageFilter.SMOOTH_MORE,
              'sharpen':ImageFilter.SHARPEN
              }

    im = toimage(arr)
    if ftype not in _tdict.keys():
        raise ValueError, "Unknown filter type."
    return fromimage(im.filter(_tdict[ftype]))


def radon(arr,theta=None):
    if theta is None:
        theta = mgrid[0:180]
    s = zeros((arr.shape[1],len(theta)), float)
    k = 0
    for th in theta:
        im = imrotate(arr,-th)
        s[:,k] = sum(im,axis=0)
        k += 1
    return s
