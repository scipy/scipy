#import Tkinter
#from Tkinter import *
#from Scientific.TkWidgets.TkPlotCanvas import *
from Numeric import *
from scipy_base.fastumath import *
#from _Graphics import *
#import PIL.ImageTk,PIL.Image
#ImageTk = PIL.ImageTk
#Image = PIL.Image
import gist, types, string, os
import Mplot
import scipy.io as io


def read_act(filename):
    fid = io.fopen('%s.act'%filename)
    p = fid.fread(256*3,'byte')
    fid.close()
    p.shape = (256,3)
    return p

Palette = {}

def make_palettes():
    p = arange(0,256,1,typecode='b')[:,NewAxis] * \
        ones((3,),typecode='b')[NewAxis,:]
    Palette['gray'] = p

    p = arange(0,256,1,typecode='b')[:,NewAxis] * \
        array([0,0,1],typecode='b')[NewAxis,:]    
    Palette['blue'] = p

    p = arange(0,256,1,typecode='b')[:,NewAxis] * \
        array([1,0,0],typecode='b')[NewAxis,:]        
    Palette['red'] = p

    p = arange(0,256,1,typecode='b')[:,NewAxis] * \
        array([0,1,0],typecode='b')[NewAxis,:]        
    Palette['green'] = p

    p = zeros((256,3),'b')
    dp = (256-40)/128.0
    p[:128,2] = arange(40,256,dp,typecode='b')[::-1]
    p[128:,0] = arange(40,256,dp,typecode='b')
    Palette['wave'] = p

    Palette['Aaasmo'] = read_act('%s/Research/Library/Python/Graphics/Aaasmo'% os.environ['HOME'])
    Palette['awave'] = read_act('%s/Research/Library/Python/Graphics/awave'% os.environ['HOME'])

def _interpolate_colormap(p,num):
    assert(len(p.shape)==2)
    assert(p.shape[1]<=3)
    pf = p.astype('d')
    pfnew = signaltools.resample(pf,num,axis=0)
    return pfnew.astype('b')
    
def read_gist_palettes():
    import commands


make_palettes()
read_gist_palettes()

##def NumPyToTkImage(data,size):
##    data = ZeroOrderHold(data,size)
##    shape = tuple(array(data.shape)[::-1])
##    return ImageTk.PhotoImage(Image.fromstring("L",shape,data.tostring()))

def _makebytes(data,top,cmin,cmax):
    if data.typecode == UnsignedInt8:
        return data
    if cmin is None:
        if cmax is None:
            bytedata = gist.bytscl(data,top=255)
        else:
            bytedata = gist.bytscl(data,cmax=cmax,top=255)
    else:
        if cmax is None:
            bytedata = gist.bytscl(data,cmin=cmin,top=255)
        else:
            bytedata = gist.bytscl(data,cmin=cmin,cmax=cmax,top=255)
    return bytedata

def _makewavebytes(data,bot,top,cmin,cmax):
    if top is None:
        top = 255
    if bot is None:
        bot = 1
    mid = (top+1) / 2        # Zero will get mapped to here.
    if cmin is None:
        cmin = min(ravel(data))
    if cmax is None:
        cmax = max(ravel(data))

    bytelow = (mid-bot)/(0.0-cmin) * data + mid
    bytehigh = (top - mid)/(cmax-0.0) * data + mid
    bytedata = where(greater(data,0),bytehigh,bytelow).astype('b')
    
    return bytedata

    
##def array2image(data,expand=None,top=255,cmin=None,cmax=None,p=None):
##    """array2image(data,expand=None,cmin=None,cmax=None,p=None)"""
##    assert(len(data.shape) == 2)
##    shape = (data.shape[1],data.shape[0])
##    bytedata = _makebytes(data,top,cmin,cmax)
##    if expand is not None:
##        if type(expand) is types.IntType:
##            expand = (expand,expand)
##        assert(type(expand) in [types.TupleType, ArrayType, types.ListType])
##        assert(expand[0] > 0 and expand[1] > 0)
##        size = tuple(array(shape)*array(expand))
##        bytedata = ZeroOrderHold(bytedata,size)
##        shape = (size[1],size[0])
##    image = Image.fromstring("L",shape,bytedata.tostring())
##    if p is not None:
##        image.putpalette(asarray(p,typecode='b').tostring())
##    return image

def write_frames(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None):
    """filelist =
    write_frames(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None)
    """
    assert(dim < 3)
    assert(len(data.shape) == 3)
    slobj = [slice(None)]*3
    filelist = []
    #data = _makebytes(data,255,cmin,cmax)  # scale all frames the same
    for k in range(data.shape[dim]):
        slobj[dim] = k
        im = array2image(data[slobj],expand,cmin=cmin,cmax=cmax,p=p)
        file = '%s%03d.png'%(filename,k)
        im.save(file)
        filelist.append(file)
    return filelist


def make_movie(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None):
    """make_movie(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None)"""
    assert(dim < 3)
    assert(len(data.shape) == 3)
    slobj = [slice(None)]*3
    filelist = []
    data = _makebytes(data,255,cmin,cmax)  # scale all frames the same
    for k in range(data.shape[dim]):
        slobj[dim] = k
        im = array2image(data[slobj],expand,top=255,cmin=cmin,cmax=cmax,p=p)
        im2 = im.convert('RGB')
        file = '/tmp/%s%03d.tiff'%(filename,k)
        im2.save(file)
        filelist.append(file)
    os.system('convert %s %s.mpg' % (string.join(filelist," "),filename))
    os.system('rm %s' % string.join(filelist," "))


def make_wavemovie(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None):
    """make_wavemovie(data,filename,dim=0,expand=None,cmin=None,cmax=None,p=None)"""
    assert(dim < 3)
    assert(len(data.shape) == 3)
    slobj = [slice(None)]*3
    filelist = []
    data = _makewavebytes(data,0,255,cmin,cmax)  # scale all frames the same
    for k in range(data.shape[dim]):
        slobj[dim] = k
        im = array2image(data[slobj],expand,top=255,cmin=cmin,cmax=cmax,p=p)
        file = '/tmp/%s%03d.gif'%(filename,k)
        im.save(file)
        filelist.append(file)
    os.system('convert %s %s.gif' % (string.join(filelist," "),filename))
    os.system('rm %s' % string.join(filelist," "))












