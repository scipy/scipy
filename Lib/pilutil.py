# Functions which need the PIL

from scipy_base import ppimport
import types
import Numeric

from scipy_base import exp, amin, amax, ravel, asarray, cast, arange, \
     ones, NewAxis, transpose, mgrid, iscomplexobj, sum, zeros

Image = ppimport('Image')
ImageFilter = ppimport('ImageFilter')

__all__ = ['fromimage','toimage','imsave','imread','bytescale',
           'imrotate','imresize','imshow','imfilter','radon']

_UInt8 = Numeric.UnsignedInt8

# Returns a byte-scaled image
def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    if data.typecode == _UInt8:
        return data
    high = high - low
    if cmin is None:
        cmin = amin(ravel(data))
    if cmax is None:
        cmax = amax(ravel(data))
    scale = high *1.0 / (cmax-cmin or 1)
    bytedata = ((data*1.0-cmin)*scale + 0.4999).astype(_UInt8)
    return bytedata + cast[_UInt8](low)
            
def imread(name,flatten=0):
    """Read an image file from a filename.

    Optional arguments:

     - flatten (0): if true, the image is flattened by calling convert('F') on
     the resulting image object.  This flattens the color layers into a single
     grayscale layer.
    """

    im = Image.open(name)
    return fromimage(im,flatten=flatten)

def imsave(name, arr):
    """Save an array to an image file.
    """
    im = toimage(arr)
    im.save(name)
    return

def fromimage(im, flatten=0):
    """Takes a PIL image and returns a copy of the image in a Numeric container.
    If the image is RGB returns a 3-dimensional array:  arr[:,:,n] is each channel

    Optional arguments:

    - flatten (0): if true, the image is flattened by calling convert('F') on
    the image object before extracting the numerical data.  This flattens the
    color layers into a single grayscale layer.  Note that the supplied image
    object is NOT modified.
    """
    assert Image.isImageType(im), "Not a PIL image."
    if flatten:
        im = im.convert('F')
    mode = im.mode
    adjust = 0
    if mode == '1':
        im = im.convert(mode='L')
        mode = 'L'
        adjust = 1
    str = im.tostring()
    type = 'b'
    if mode == 'F':
        type = 'f'
    if mode == 'I':
        type = 'i'
    arr = Numeric.fromstring(str,type)
    shape = list(im.size)
    shape.reverse()
    if mode == 'P':
        arr.shape = shape
        if im.palette.rawmode != 'RGB':
            print "Warning: Image has invalid palette."
            return arr
        pal = Numeric.fromstring(im.palette.data,type)
        N = len(pal)
        pal.shape = (int(N/3.0),3)
        return arr, pal
    if mode in ['RGB','YCbCr']:
        shape += [3]
    elif mode in ['CMYK','RGBA']:
        shape += [4]
    arr.shape = shape
    if adjust:
        arr = (arr != 0)
    return arr

_errstr = "Mode is unknown or incompatible with input array shape."
def toimage(arr,high=255,low=0,cmin=None,cmax=None,pal=None,
            mode=None,channel_axis=None):
    """Takes a Numeric array and returns a PIL image.  The mode of the
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

    The Numeric array must be either 2 dimensional or 3 dimensional.
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
            image = Image.fromstring(mode,shape,data.astype('f').tostring())
            return image
        if mode in [None, 'L', 'P']:
            bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
            image = Image.fromstring('L',shape,bytedata.tostring())
            if pal is not None:
                image.putpalette(asarray(pal,typecode=_UInt8).tostring())
                # Becomes a mode='P' automagically.
            elif mode == 'P':  # default gray-scale
                pal = arange(0,256,1,typecode='b')[:,NewAxis] * \
                      ones((3,),typecode='b')[NewAxis,:]
                image.putpalette(asarray(pal,typecode=_UInt8).tostring())
            return image
        if mode == '1':  # high input gives threshold for 1
            bytedata = ((data > high)*255).astype('b')
            image = Image.fromstring('L',shape,bytedata.tostring())   
            image = image.convert(mode='1')
            return image
        if cmin is None:
            cmin = amin(ravel(data))
        if cmax is None:
            cmax = amax(ravel(data))
        data = (data*1.0 - cmin)*(high-low)/(cmax-cmin) + low
        if mode == 'I':
            image = Image.fromstring(mode,shape,data.astype('i').tostring())
        else:
            raise ValueError, _errstr
        return image

    # if here then 3-d array with a 3 or a 4 in the shape length.
    # Check for 3 in datacube shape --- 'RGB' or 'YCbCr'
    if channel_axis is None:
        if (3 in shape):
            ca = Numeric.nonzero(asarray(shape) == 3)[0]
        else:
            ca = Numeric.nonzero(asarray(shape) == 4)
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
    """Rotate an image counter-clockwise by angle degrees.

    Interpolation methods can be:
        'nearest' :  for nearest neighbor
        'bilinear' : for bilinear
        'cubic' or 'bicubic' : for bicubic 
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
    if (len(arr.shape) == 3) and (arr.shape[2] == 4):
        try:
            import os
            im.save('/tmp/scipy_imshow.png')
            if os.system("(xv /tmp/scipy_imshow.png; rm -f /tmp/scipy_imshow.png)&"):
                raise RuntimeError
            return
        except:
            print "Warning: Alpha channel may not be handled correctly."
            
    im.show()
    return

def imresize(arr,size):
    """Resize an image.

    If size is an integer it is a percentage of current size.
    If size is a float it is a fraction of current size.
    If size is a tuple it is the size of the output image.
    """
    im = toimage(arr)
    ts = type(size)
    if ts is types.IntType:
        size = size / 100.0
    if type(size) is types.FloatType:
        size = (im.size[0]*size,im.size[1]*size)
    else:
        size = (size[1],size[0])
    imnew = im.resize(size)
    return fromimage(imnew)


def imfilter(arr,ftype):
    """Simple filtering of an image.

    type can be:
            'blur', 'contour', 'detail', 'edge_enhance', 'edge_enhance_more',
            'emboss', 'find_edges', 'smooth', 'smooth_more', 'sharpen'
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
    s = zeros((arr.shape[1],len(theta)),'d')
    k = 0
    for th in theta:
        im = imrotate(arr,-th)
        s[:,k] = sum(im,axis=0)
        k += 1
    return s
        
