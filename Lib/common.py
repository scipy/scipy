# Functions which are common and require SciPy Base and Level 1 SciPy
# (stats, special, linalg)

# Needs eigenvalues
import Numeric
import types, sys
from scipy import special, stats
from scipy_base import exp, amin, amax, ravel, asarray, cast, arange, \
     ones, NewAxis, transpose
import scipy_base.fastumath

__all__ = ['factorial','comb','rand','randn','disp','who','bytescale']
    
def factorial(n,exact=0):
    """n! = special.gamma(n+1)

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed."""
    if n < 0:
        raise ValueError, "n must be >= 0"
    if exact:
        n = int(n)
        val = 1L
        for k in xrange(1,n+1):
            val = val*k
        return val
    else:
        return special.gamma(n+1)


def comb(N,k,exact=0):
    """Combinations of N things taken k at a time.

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed.
    """
    if (k > N) or (N < 0) or (k < 0):
        raise ValueError, "N and k must be non-negative and k <= N"
    if exact:
        N,k = map(int,(N,k))
        val = 1L
        for n in xrange(N-k+1,N+1):
            val = val*n
        for n in xrange(1,k+1):
            val = val / n
        return val
    else:
        lgam = special.gammaln
        return exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))


def rand(*args):
    """rand(d1,...,dn) returns a matrix of the given dimensions
    which is initialized to random numbers from a uniform distribution
    in the range [0,1).
    """
    return stats.random(args)

def randn(*args):
    """u = randn(d0,d1,...,dn) returns zero-mean, unit-variance Gaussian
    random numbers in an array of size (d0,d1,...,dn)."""
    return stats.stnorm(size=args)

def lena():
    import cPickle, os
    d,junk = os.path.split(os.path.abspath(scipy.__file__))
    fname = os.path.join(d,'plt','lena.dat')
    f = open(fname,'rb')
    lena = scipy.array(cPickle.load(f))
    f.close()
    return lena

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
    scale = high *1.0 / (cmax-cmin)
    bytedata = ((data*1.0-cmin)*scale + 0.4999).astype(_UInt8)
    return bytedata + cast[_UInt8](low)
            
try:
    import Image
    def fromimage(im):
        """Takes a PIL image and returns a copy of the image in a Numeric container.
        If the image is RGB returns a 3-dimensional array:  arr[:,:,n] is each channel
        """
        assert Image.isImageType(im), "Not a PIL image."
        mode = im.mode
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
            shape = (3,) + shape
        elif mode in ['CMYK','RGBA']:
            shape = (4,) + shape
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
        data = Numeric.asarray(arr)
        shape = list(data.shape)
        valid = len(shape)==2 or ((len(shape)==3) and \
                                  ((3 in shape) or (4 in shape)))
        assert valid, "Not a suitable array shape for any mode."
        if len(shape) == 2:
            shape = (shape[1],shape[0]) # columns show up first
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
            if mode == 'F':
                image = Image.fromstring(mode,shape,data.astype('f').tostring())
            elif mode == 'I':
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
            
        if shape[ca] not in [3,4]:
            raise ValueError, "Channel axis dimension is not valid."

        numch = shape[ca]
        bytedata = bytescale(data,high=high,low=low,cmin=cmin,cmax=cmax)
        if ca == 0:
            strdata = bytedata.tostring()
            shape = (shape[2],shape[1])
        elif ca == 1:
            strdata = transpose(bytedata,(1,0,2)).tostring()
            shape = (shape[2],shape[0])
        elif ca == 2:
            strdata = transpose(bytedata,(2,0,1)).tostring()
            shape = (shape[1],shape[0])
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

    __all__.extend(['fromimage','toimage'])

except ImportError:
    pass
            
                                  
                                 
        
    

#-----------------------------------------------------------------------------
# Matlab like functions for output and information on the variables used.
#-----------------------------------------------------------------------------
def disp(mesg, device=None, linefeed=1):
    """Display a message to device (default is sys.stdout) with(out) linefeed.
    """
    if device is None:
        device = sys.stdout
    if linefeed:
        device.write('%s\n' % mesg)
    else:
        device.write('%s' % mesg)
    device.flush()
    return

def who(vardict=None):
    """Print the Numeric arrays in the given dictionary (or globals() if None).
    """
    if vardict is None:
        print "Pass in a dictionary:  who(globals())"
        return
    sta = []
    cache = {}
    for name in vardict.keys():
        if isinstance(vardict[name],Numeric.ArrayType):
            var = vardict[name]
            idv = id(var)
            if idv in cache.keys():
                namestr = name + " (%s)" % cache[idv]
                original=0
            else:
                cache[idv] = name
                namestr = name
                original=1
            shapestr = " x ".join(map(str, var.shape))
            bytestr = str(var.itemsize()*Numeric.product(var.shape))
            sta.append([namestr, shapestr, bytestr, _namefromtype[var.typecode()], original])

    maxname = 0
    maxshape = 0
    maxbyte = 0
    totalbytes = 0
    for k in range(len(sta)):
        val = sta[k]
        if maxname < len(val[0]):
            maxname = len(val[0])
        if maxshape < len(val[1]):
            maxshape = len(val[1])
        if maxbyte < len(val[2]):
            maxbyte = len(val[2])
        if val[4]:
            totalbytes += int(val[2])

    max = Numeric.maximum
    if len(sta) > 0:
        sp1 = max(10,maxname)
        sp2 = max(10,maxshape)
        sp3 = max(10,maxbyte)
        prval = "Name %s Shape %s Bytes %s Type" % (sp1*' ', sp2*' ', sp3*' ')
        print prval + "\n" + "="*(len(prval)+5) + "\n"
        
    for k in range(len(sta)):
        val = sta[k]
        print "%s %s %s %s %s %s %s" % (val[0], ' '*(sp1-len(val[0])+4),
                                        val[1], ' '*(sp2-len(val[1])+5),
                                        val[2], ' '*(sp3-len(val[2])+5),
                                        val[3])
    print "\nUpper bound on total bytes  =       %d" % totalbytes
    return
    

#-----------------------------------------------------------------------------

def test(level=10):
    from scipy_base.testing import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_base.testing import module_test_suite
    return module_test_suite(__name__,__file__,level=level)

