# 1998 - 2003
# Author: Travis Oliphant
# Copyright: SciPy

import gist
import pl3d, plwf
import Numeric
from Numeric import ravel, reshape, repeat, arange, transpose, compress, \
     where, ones, NewAxis, asarray
import MLab
from MLab import pi, cos, sin, arctan2, array, angle
import types
import write_style
from write_style import inches, points
import scipy
import scipy_base
import scipy.signal as signal

_dpi = 75
_hold = 0
_maxcolors=256
_textcolor=None

gist.pldefault(dpi=_dpi,maxcolors=_maxcolors)

# Get a directory that the user has write access to for
#  storing temporary *.gp and *.gs files
import tempfile
import os
import sys

# functions copied from weave.catalog

def whoami():
    """return a string identifying the user."""
    return os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"

def create_dir(p):
    """ Create a directory and any necessary intermediate directories."""
    if not os.path.exists(p):
        try:
            os.mkdir(p)
        except OSError:
            # perhaps one or more intermediate path components don't exist
            # try to create them
            base,dir = os.path.split(p)
            create_dir(base)
            # don't enclose this one in try/except - we want the user to
            # get failure info
            os.mkdir(p)

def is_writable(dir):
    dummy = os.path.join(dir, "dummy")
    try:
        open(dummy, 'w')
    except IOError:
        return 0
    os.unlink(dummy)
    return 1

# end functions copied from weave.catalog

def _getdir(name='scipy_xplt'):
    try:
        path = os.path.join(os.environ['HOME'],'.' + name)
    except KeyError:
        path = os.path.join(tempfile.gettempdir(),"%s"%whoami(),
                            name)
    if not os.path.exists(path):
        create_dir(path)
        os.chmod(path,0700)
    if not is_writable(path):
        print "warning: default directory is not write accessible."
        print "default:", path
    return path

_user_path = _getdir()

try:
    import Scientific.Statistics.Histogram
    SSH = Scientific.Statistics.Histogram
    def histogram(data,nbins=80,range=None,ntype=0,bar=1,bwidth=0.8,bcolor=0):
        """Plot a histogram.  ntype is the normalization type.

        Use ntype == 2 to compare with probability density function.
        """
        h = SSH.Histogram(data,nbins,range)
        if ntype == 1:
            h.normalize()
        elif ntype == 2:
            h.normalizeArea()
        if bar:
            barplot(h[:,0],h[:,1],width=bwidth,color=bcolor)
        else:
            plot(h[:,0],h[:,1])
        return h
except ImportError:
    try:
        import Statistics
        SSH = Statistics
        def histogram(data,nbins=80,range=None,ntype=0,bar=1,bwidth=0.8,bcolor=0):
            """Plot a histogram.  ntype is the normalization type.
            
            Use ntype == 2 to compare with probability density function.
            """
            h = SSH.histogram(data,nbins,range)
            if ntype == 1:
                h.normalize()
            elif ntype == 2:
                h.normalizeArea()
            if bar:
                barplot(h[:,0],h[:,1],width=bwidth,color=bcolor)                
            else:
                plot(h[:,0],h[:,1])
            return h        
    except ImportError:
        from scipy.stats import histogram2 as hist
        def histogram(data,nbins=80,range=None,ntype=0,bar=1,bwidth=0.8,bcolor=0):
            """Plot a histogram.  ntype is the normalization type.
            
            Use ntype == 2 to compare with probability density function.
            """
            if range is None:
                dmin = Numeric.minimum.reduce(data)
                dmax = Numeric.maximum.reduce(data)
            else:
                dmin, dmax = range
            dmin = dmin + 0.0
            dmax = dmax + 0.0
            bin_width = (dmax - dmin)/nbins
            darray = Numeric.zeros((nbins,2),Numeric.Float)
            darray[:,0] = dmin + bin_width*(Numeric.arange(nbins)+0.5)
            bins = dmin + bin_width*(Numeric.arange(nbins))
            darray[:,1] = hist(data,bins)
            if ntype == 1:
                darray[:,1] = 1.0*darray[:,1] / Numeric.add.reduce(darray[:,1])
            elif ntype == 2:
                darray[:,1] = 1.0/bin_width*darray[:,1] / \
                              Numeric.add.reduce(darray[:,1])
            if bar:
                barplot(darray[:,0],darray[:,1],width=bwidth,color=bcolor)
            else:
                plot(darray[:,0],darray[:,1])
            return darray


def textcolor(color=None):
    global _textcolor
    if color is not None:
        _textcolor = color
    return _textcolor
    
def reverse_dict(dict):
    newdict = {}
    for key in dict.keys():
        newdict[dict[key]] = key
    return newdict

_colornum = {'black':-3, 'white':-4,'red':-5,'green':-6,'blue':-7,'cyan':-8,'magenta':-9,'yellow':-10}
_types = {'-':'solid','|':'dash',':':'dot','-.':'dashdot','-:':'dashdotdot'}
_corder = ['b','r','m','g','c','k','y']
_colors = {'k':'black','r':'red','b':'blue','m':'magenta','g':'green','y':'yellow','c':'cyan','w':'white'}
_markers = { '+':'\2','.':'\1','*':'\3','o':'\4','x':'\5'}
_current_style='work.gs'

_rtypes = reverse_dict(_types)
_rtypes['none'] = ''
_rcolors = reverse_dict(_colors)
_rmarkers = reverse_dict(_markers)

def _find_and_set(dict, str, default):
    import string
    value = default
    for k in dict.keys():
        if string.find(str,k) >= 0:
            value = dict[k]
            break
    return value

def barplot(x,y,width=0.8,color=0):
    """Plot a barplot.

    Description:
    
      Plot a barplot with centers at x and heights y with given color

    Inputs:

      x, y -- Centers and heights of bars
      width -- Relative width of the bars.
      color -- A number from the current palette.
    """
    N = 4*Numeric.ones(len(x))
    hw = width * (x[1]-x[0])/ 2.0
    Xa = x-hw
    Xb = x+hw
    Ya = Numeric.zeros(len(y),'d')
    Yb = y
    X = Numeric.array((Xa,Xa,Xb,Xb))
    Y = Numeric.array((Ya,Yb,Yb,Ya))
    X = Numeric.reshape(Numeric.transpose(X),(4*len(N),))
    Y = Numeric.reshape(Numeric.transpose(Y),(4*len(N),))
    try:
        override = 1
        savesys = gist.plsys(2)
        gist.plsys(savesys)
    except:
        override = 0
    if _hold or override:
        pass
    else:
        gist.fma()
    Z = color * Numeric.ones(len(N))
    gist.plfp(Z.astype('b'),Y,X,N)
    return

def hold(state):
    """Draw subsequent plots over the current plot.

    Inputs:

      state -- If 'on' or 'yes' hold the current plot.
               Otherwise refresh screen when drawing new plot.
    """
    global _hold
    if state in ['on', 'yes']:
        _hold = 1
    elif state in ['off', 'no']:
        _hold = 0
    else:
        raise ValueError, 'holds argument can be "on","off",'\
                          '"yes","no". Not ' + state
    return
          

def errorbars(x,y,err,ptcolor='r',linecolor='b',pttype='o',linetype='-',fac=0.25):
    """Draw connected points with errorbars.

    Description:

      Plot connected points with errorbars.

    Inputs:

      x, y -- The points to plot.
      err -- The error in the y values.
      ptcolor -- The color for the points.
      linecolor -- The color of the connecting lines and error bars.
      pttype -- The type of point ('o', 'x', '+', '.', 'x', '*')
      linetype -- The type of line ('-', '|', ':', '-.', '-:')
      fac -- Adjusts how long the horizontal lines are which make the
             top and bottom of the error bars.
    """
    # create line arrays
    yb = y - err
    ye = y + err
    try:
        override = 1
        savesys = gist.plsys(2)
        gist.plsys(savesys)
    except:
        override = 0
    if _hold or override:
        pass
    else:
        gist.fma()
    y = where(scipy.isfinite(y),y,0)
    gist.plg(y,x,color=_colors[ptcolor],marker=_markers[pttype],type='none')
    gist.pldj(x,yb,x,ye,color=_colors[linecolor],type=_types[linetype])
    viewp = gist.viewport()
    plotlims = gist.limits()
    conv_factorx = (viewp[1] - viewp[0]) / (plotlims[1]-plotlims[0])
    conv_factory = (viewp[3] - viewp[2]) / (plotlims[3]-plotlims[2])    
    width = fac*(x[1]-x[0])
    x0 = x-width/2.0
    x1 = x+width/2.0
    gist.pldj(x0,ye,x1,ye,color=_colors[linecolor],type=_types[linetype])
    gist.pldj(x0,yb,x1,yb,color=_colors[linecolor],type=_types[linetype])
    return

def legend(text,linetypes=None,lleft=None,color=None,tfont='helvetica',fontsize=14,nobox=0):
    """Construct and place a legend.

    Description:

      Build a legend and place it on the current plot with an interactive
      prompt.

    Inputs:

      text -- A list of strings which document the curves.
      linetypes -- If not given, then the text strings are associated
                   with the curves in the order they were originally
                   drawn.  Otherwise, associate the text strings with the
                   corresponding curve types given.  See plot for description.
                   
    """    
    global _hold
    global _textcolor
    if color is None:
        color = _textcolor
    else:
        _textcolor = color
    if color is None:
        color = 'black'

    sys = gist.plsys()
    if sys == 0:
        gist.plsys(1)
    viewp = gist.viewport()
    gist.plsys(sys)
    DX = viewp[1] - viewp[0]
    DY = viewp[3] - viewp[2]
    width = DY / 10.0;
    if lleft is None:
        lleft = gist.mouse(0,0,"Click on point for lower left coordinate.")
        llx = lleft[0]
        lly = lleft[1]
    else:
        llx,lly = lleft[:2]

    savesys = gist.plsys()
    dx = width / 3.0
    legarr = Numeric.arange(llx,llx+width,dx)
    legy = Numeric.ones(legarr.shape)
    dy = fontsize*points*1.2
    deltay = fontsize*points / 2.8
    deltax = fontsize*points / 2.6 * DX / DY
    ypos = lly + deltay;
    if linetypes is None:
        linetypes = _GLOBAL_LINE_TYPES[:]  # copy them out
    gist.plsys(0)
    savehold = _hold
    _hold = 1
    for k in range(len(text)):
        plot(legarr,ypos*legy,linetypes[k])
        print linetypes[k], text[k]
        print llx+width+deltax, ypos-deltay
        if text[k] != "":
            gist.plt(text[k],llx+width+deltax,ypos-deltay,
                     color=color,font=tfont,height=fontsize,tosys=0)
        ypos = ypos + dy
    _hold = savehold
    if nobox:
        pass
    else:
        gist.plsys(0)
        maxlen = MLab.max(map(len,text))
        c1 = (llx-deltax,lly-deltay)
        c2 = (llx + width + deltax + fontsize*points* maxlen/1.8 + deltax,
              lly + len(text)*dy)
        linesx0 = [c1[0],c1[0],c2[0],c2[0]]
        linesy0 = [c1[1],c2[1],c2[1],c1[1]]
        linesx1 = [c1[0],c2[0],c2[0],c1[0]]
        linesy1 = [c2[1],c2[1],c1[1],c1[1]]
        gist.pldj(linesx0,linesy0,linesx1,linesy1,color=color)
    gist.plsys(savesys)
    return

def ispointtype(linetype):
    if len(linetype) > 2:
        return 0
    if (len(linetype) == 1):
        if linetype[0] in _markers.keys():
            return 1
        else:
            return 0
    if linetype[0] in _markers.keys():
        if linetype[1] in _colors.keys():
            return 1
        else:
            return 0
    if linetype[0] in _colors.keys():
        if linetype[1] in _markers.keys():
            return 1
        else:
            return 0
    return 0

##def legend(text,linetypes=None,lleft=None,color='black',tfont='helvetica',fontsize=14,nobox=0):
##    viewp = gist.viewport()
##    plotlims = gist.limits()
##    gist.limits(plotlims)
##    conv_factorx = (viewp[1] - viewp[0]) / (plotlims[1]-plotlims[0])
##    conv_factory = (viewp[3] - viewp[2]) / (plotlims[3]-plotlims[2])

##    width = (plotlims[1] - plotlims[0]) / 10.0;
##    if lleft is None:
##        lleft = gist.mouse(-1,0,"Click on point for lower left coordinate.")
##        llx = lleft[0]
##        lly = lleft[1]
##    else:
##        llx,lly = lleft[:2]

##    dx = width / 3.0
##    legarr = Numeric.arange(llx,llx+width,dx)
##    legy = Numeric.ones(legarr.shape)
##    dy = fontsize*points/conv_factory*1.15
##    deltay = fontsize*points / conv_factory / 2.8
##    deltax = fontsize*points / conv_factorx / 2.8
##    ypos = lly + deltay;
##    if linetypes is None:
##        linetypes = _GLOBAL_LINE_TYPES[:]  # copy them out
##    for k in range(len(text)):
##        if ispointtype(linetypes[k]):
##            pt = len(legarr)/2
##            plot([legarr[pt]],[ypos*legy[pt]],linetypes[k], hold=1)
##        else:
##            plot(legarr,ypos*legy,linetypes[k],hold=1)
##        print llx+width+deltax, ypos-deltay
##        if text[k] != "":
##            gist.plt(text[k],llx+width+deltax,ypos-deltay,
##                     color=color,font=tfont,height=fontsize,tosys=1)
##        ypos = ypos + dy

##    if nobox:
##        pass
##    else:
##        maxlen = MLab.max(map(len,text))
##        c1 = (llx-deltax,lly-deltay)
##        c2 = (llx + width + deltax + fontsize*points/conv_factorx * maxlen/1.8 + deltax,
##              lly + len(text)*dy)
##        linesx0 = [c1[0],c1[0],c2[0],c2[0]]
##        linesy0 = [c1[1],c2[1],c2[1],c1[1]]
##        linesx1 = [c1[0],c2[0],c2[0],c1[0]]
##        linesy1 = [c2[1],c2[1],c1[1],c1[1]]
##        gist.pldj(linesx0,linesy0,linesx1,linesy1,color=color)
##    return

import operator
def arrow(x0,y0,x1,y1,color=0,ang=45.0,height=6,width=1.5,lc=None):
    """Draw an arrow.   

    Description:

      Draw an arrow from (x0,y0) to (x1,y1) in the current coordinate system.

    Inputs:

      x0, y0 -- The beginning point.
      x1, y1 -- Then ending point.
      color -- The color of the arrowhead.  Number represents an index
               in the current palette or a negative number or a spelled
               out basic color. 
      lc -- The color of the line (same as color by default).
      ang -- The angle of the arrowhead.
      height -- The height of the arrowhead in points.
      width -- The width of the arrow line in points.
    """
    if lc is None:
        lc = color
    if type(lc) is types.StringType:
        lc = _colornum[lc]
    if type(color) is types.StringType:
        color = _colornum[color]
    vp = gist.viewport()
    plotlims = gist.limits()
    gist.limits(plotlims)
    conv_factorx = (vp[1]-vp[0]) / (plotlims[1]-plotlims[0])
    conv_factory = (vp[3]-vp[2]) / (plotlims[3]-plotlims[2])
    ang = ang*pi/180
    height = height*points
    hypot = height / cos(ang)
    difx = (x1 - x0) * conv_factorx
    dify = (y1 - y0) * conv_factory
    theta = arctan2(dify,difx) + pi
    tha = theta + ang
    thb = theta - ang
    x1a = x1 + hypot*cos(tha) / conv_factorx
    x1b = x1 + hypot*cos(thb) / conv_factorx
    y1a = y1 + hypot*sin(tha) / conv_factory
    y1b = y1 + hypot*sin(thb) / conv_factory
    gist.pldj([x0],[y0],[x1],[y1],color=lc,width=width)
    gist.plfp(array([color],'b'),[y1,y1a,y1b],[x1,x1a,x1b],[3])
    return

def _parse_type_arg(thearg,nowplotting):
    indx = nowplotting % len(_corder)
    if type(thearg) is type(''):
        tomark = 1

        thetype = _find_and_set(_types,thearg,'none')
        thecolor = _find_and_set(_colors,thearg,_colors[_corder[indx]])
        themarker = _find_and_set(_markers,thearg,None)
        
        if (themarker is None):
            tomark = 0
            if thetype == 'none':
                thetype = 'solid'        

        return (thetype, thecolor, themarker, tomark)

    else:  # no string this time
        return ('solid',_colors[_corder[indx]],'Z',0)

_GLOBAL_LINE_TYPES=[]
def clear_global_linetype():
    for k in range(len(_GLOBAL_LINE_TYPES)):
        _GLOBAL_LINE_TYPES.pop()

def append_global_linetype(arg):
    _GLOBAL_LINE_TYPES.append(arg)

def _minsqueeze(arr,min=1):
    # eliminate extra dimensions above min
    arr = asarray(arr)
    arr = scipy.squeeze(arr)
    n = len(arr.shape)
    if n < min:
        arr.shape = arr.shape + (1,)*(min-n)
    return arr

def plot(x,*args,**keywds):
    """Plot curves.

    Description:

      Plot one or more curves on the same graph.

    Inputs:

      There can be a variable number of inputs which consist of pairs or
      triples.  The second variable is plotted against the first using the
      linetype specified by the optional third variable in the triple.  If
      only two plots are being compared, the x-axis does not have to be
      repeated.
    """
    try:
        override = 1
        savesys = gist.plsys(2)
        gist.plsys(savesys)
    except:
        override = 0
    global _hold
    try: _hold=keywds['hold']
    except KeyError: pass
    try: linewidth=float(keywds['width'])
    except KeyError: linewidth=1.0
    try: msize = float(keywds['msize'])
    except KeyError: msize=1.0
    if _hold or override:
        pass
    else:
        gist.fma()
    gist.animate(0)
    savesys = gist.plsys()
    winnum = gist.window()
    if winnum < 0:
        gist.window(0)
    if savesys >= 0:
        gist.plsys(savesys)
    nargs = len(args)
    if nargs == 0:
        y = _minsqueeze(x)
        x = Numeric.arange(0,len(y))
        if scipy.iscomplexobj(y):
            print "Warning: complex data plotting real part."
            y = y.real
        y = where(scipy.isfinite(y),y,0)
        gist.plg(y,x,type='solid',color='blue',marks=0,width=linewidth)
        return
    y = args[0]
    argpos = 1
    nowplotting = 0
    clear_global_linetype()
    while 1:
        try:
            thearg = args[argpos]
        except IndexError:
            thearg = 0
        thetype,thecolor,themarker,tomark = _parse_type_arg(thearg,nowplotting)
        if themarker == 'Z':  # args[argpos] was data or non-existent.
            pass
            append_global_linetype(_rtypes[thetype]+_rcolors[thecolor])
        else:                 # args[argpos] was a string
            argpos = argpos + 1
            if tomark:
                append_global_linetype(_rtypes[thetype]+_rcolors[thecolor]+_rmarkers[themarker])
            else:
                append_global_linetype(_rtypes[thetype]+_rcolors[thecolor])
        if scipy.iscomplexobj(x) or scipy.iscomplexobj(y):
            print "Warning: complex data provided, using only real part."
            x = scipy.real(x)
            y = scipy.real(y)
	y = where(scipy.isfinite(y),y,0)
        y = _minsqueeze(y) 
        x = _minsqueeze(x)
        gist.plg(y,x,type=thetype,color=thecolor,marker=themarker,marks=tomark,msize=msize,width=linewidth)

        nowplotting = nowplotting + 1

        ## Argpos is pointing to the next potential triple of data.
        ## Now one of four things can happen:
        ##
        ##   1:  argpos points to data, argpos+1 is a string
        ##   2:  argpos points to data, end
        ##   3:  argpos points to data, argpos+1 is data
        ##   4:  argpos points to data, argpos+1 is data, argpos+2 is a string

        if argpos >= nargs: break      # no more data

        if argpos == nargs-1:          # this is a single data value.
            x = x
            y = args[argpos]
            argpos = argpos+1
        elif type(args[argpos+1]) is types.StringType:
            x = x
            y = args[argpos]
            argpos = argpos+1
        else:   # 3 
            x = args[argpos]
            y = args[argpos+1]
            argpos = argpos+2
    return

def matplot(x,y=None,axis=-1):
    if y is None:   # no axis data
        y = x
        x = Numeric.arange(0,y.shape[axis])
    x,y = Numeric.asarray(x), Numeric.asarray(y)
    assert(len(y.shape)==2)
    assert(len(x)==y.shape[axis])
    otheraxis = (1+axis) % 2
    sliceobj = [slice(None)]*2
    if not _hold and gist.plsys() < 2:
        gist.fma()
    clear_global_linetype()
    for k in range(y.shape[otheraxis]):
        thiscolor = _colors[_corder[k % len(_corder)]] 
        sliceobj[otheraxis] = k
	ysl = where(scipy.isfinite(y[sliceobj]),y[sliceobj],0)
        gist.plg(ysl,x,type='solid',color=thiscolor,marks=0)
        append_global_linetype(_rcolors[thiscolor]+'-')


def addbox(x0,y0,x1,y1,color='black',width=1,type='-'):
    if not isinstance(color,types.IntType):
        color = _colornum[color]
    wordtype = _types[type]
    gist.pldj([x0,x1,x1,x0],[y0,y0,y1,y1],[x1,x1,x0,x0],[y0,y1,y1,y0],
              color=color,type=wordtype,width=width)

def write_palette(tofile,pal):
    pal = Numeric.asarray(pal)
    if pal.typecode() not in ['b','1','s','i','l']:
        raise ValueError, "Palette data must be integer data."
    palsize = pal.shape
    if len(palsize) > 2:
        raise TypeError, "Input must be a 1-d or 2-d array"
    if len(palsize) == 2:
        if palsize[0] == 1 and palsize[1] > 1:
            pal = pal[0]
        if palsize[1] == 1 and palsize[0] > 1:
            pal = pal[:,0]
        palsize = pal.shape
    if len(palsize) == 1:
        pal = Numeric.multiply.outer(pal,ones((3,),pal.typecode()))
        palsize = pal.shape
    if not (palsize[1] == 3 or palsize[0] == 3):
        raise TypeError, "If input is 2-d, the length of at least one dimension must be 3."

    if palsize[0] == 3 and palsize[1] != 3:
        pal = Numeric.transpose(pal)
        palsize = pal.shape

    if palsize[0] > 256:
        raise ValueError, "Palettes should be no longer than 256."
    fid = open(tofile,'w')
    fid.write("ncolors=%d\n\n#  r   g   b\n" % palsize[0])
    for k in range(palsize[0]):
        fid.write("%4d%4d%4d\n" % tuple(pal[k]))
    fid.close()

def list_palettes():
    import os, glob
    direc = os.environ['GISTPATH']
    files = glob.glob1(direc,"*.gp")
    lengths = map(len,files)
    maxlen = scipy.amax(lengths)
    print "Available palettes..."
    print "=====================\n"
    for file in files:
        print file[:-3] + ' '*(maxlen-len(file[:-3])-3) + ' --- ',
        k = 0
        fid = open(direc+"/"+file)        
        while 1:
            line = fid.readline()
            if line[0] != '#':
                fid.close()
                if k == 0:
                    print
                break
            if k > 0:
                print ' '*(maxlen+3) + line[1:-1]
            else:
                print line[1:-1]
            k = k + 1

                                     
def change_palette(pal):
    if pal is not None:
        if isinstance(pal, types.StringType):
            try:
                gist.palette('%s.gp' % pal)
            except gist.error:
                if len(pal) > 3 and pal[-2:] == 'gp':
                    gist.palette(pal)
                else:
                    raise ValueError, "Palette %s not found." % pal
        else:
            data = Numeric.transpose(Numeric.asarray(pal))
	    data = data.astype('b')
	    gist.palette(*transpose(data))
            #filename = os.path.join(_user_path,'_temp.gp')
            #write_palette(filename,data)
            #gist.palette(filename)

chpal = change_palette

def matview(A,cmax=None,cmin=None,palette=None,color='black'):
    """Plot an image of a matrix.
    """
    A = Numeric.asarray(A)
    if A.typecode() in ['D','F']:
        print "Warning: complex array given, plotting magnitude."
        A = abs(A)
    M,N = A.shape
    A = A[::-1,:]
    if cmax is None:
        cmax = max(ravel(A))
    if cmin is None:
        cmin = min(ravel(A))
    cmax = float(cmax)
    cmin = float(cmin)
    byteimage = gist.bytscl(A,cmin=cmin,cmax=cmax)
    change_palette(palette)
    gist.window(style='nobox.gs')
    _current_style='nobox.gs'
    gist.pli(byteimage)
    old_vals = gist.limits(square=1)
    vals = gist.limits(square=1)
    if color is None:
        return
    vp = gist.viewport()
    axv,bxv,ayv,byv = vp
    axs,bxs,ays,bys = vals[:4]
    # bottom left corner column
    posy = -ays*(byv-ayv)/(bys-ays) + ayv
    posx = -axs*(bxv-axv)/(bxs-axs) + axv
    gist.plt('1',posx,posy-0.005,justify='LT',color=color)
    # bottom left corner row
    gist.plt(str(M),posx-0.005,posy,justify='RB',color=color)
    # top left corner row
    posy = (M-ays)*(byv-ayv)/(bys-ays) + ayv
    gist.plt('1',posx-0.005,posy,justify='RT',color=color)
    # bottom right column
    posy = -ays*(byv-ayv)/(bys-ays) + ayv
    posx = (N-axs)*(bxv-axv)/(bxs-axs) + axv
    gist.plt(str(N),posx,posy-0.005,justify='RT',color=color)


def imagesc(z,cmin=None,cmax=None,xryr=None,_style='default', palette=None,
            color='black',colormap=None):
    """Plot an image on axes.

    z -- The data
    cmin -- Value to map to lowest color in palette (min(z) if None)
    cmax -- Value to map to highest color in palette (max(z) if None)
    xryr -- (xmin, ymin, xmax, ymax) coordinates to print
            (0, 0, z.shape[1], z.shape[0]) if None
    _style -- A 'style-sheet' to use if desired (a default one will be used
              if 'default').  If None, then no style will be imposed.
    palette -- A string for a palette previously saved in a file (see write_palette)
               or an array specifying the red-green-blue values (2-d array N x 3) or
               gray-scale values (2-d array N x 1 or 1-d array).
    color -- The color to use for the axes.
    """
    if xryr is None:
        xryr = (0,0,z.shape[1],z.shape[0])
    try:
        _style = None
        saveval = gist.plsys(2)
        gist.plsys(saveval)
    except:
        _style = 'default'        
        if not _hold:
            gist.fma()
        gist.animate(0)
    
    if _style is not None:
        if _style == "default":
            _style=os.path.join(_user_path,'image.gs')
            system = write_style.getsys(hticpos='below',vticpos='left',frame=1,
                                        color=color)
            fid = open(_style,'w')
            fid.write(write_style.style2string(system))
            fid.close()
        gist.window(style=_style)
        _current_style=_style
    if cmax is None:
        cmax = max(ravel(z))
    if cmin is None:
        cmin = min(ravel(z))
    cmax = float(cmax)
    cmin = float(cmin)
    byteimage = gist.bytscl(z,cmin=cmin,cmax=cmax)
    if (colormap is not None): palette=colormap
    change_palette(palette)
    gist.pli(byteimage,xryr[0],xryr[1],xryr[2],xryr[3])
    return


def movie(data,aslice,plen,loop=1,direc='z',cmax=None,cmin=None):
    "movie(data,slice,pause,loop=1,direc='z')"
    gist.animate(1)
    if type(aslice) is types.IntType:
        num = aslice
        aslice = [slice(None)]*3
        aslice[ord('x')-ord(direc)-1] = num
    for num in range(loop):
        for k in range(data.shape[0]):
            gist.fma()
            gist.pli(data[k][aslice],cmax=cmax,cmin=cmin)
            gist.pause(plen)
    gist.animate(0)

def setdpi(num):
    """ Set the dpi for new windows """
    if num in [75,100]:
        _dpi = num
        gist.set_default_dpi(_dpi)
    else:
        raise ValueError, "DPI must be 75 or 100"

def figure(n=None,style=os.path.join(_user_path,"currstyle.gs"), color=-2, frame=0, labelsize=14, labelfont='helvetica',aspect=1.618):
    global _figures
    if (aspect < 0.1) or (aspect > 10):
        aspect = 1.618
    if isinstance(color, types.StringType):
        color = _colornum[color]
    fid = open(style,'w')
    syst = write_style.getsys(color=color,frame=frame,
                             labelsize=labelsize,font=labelfont)
    cntr = (5.5*inches,4.25*inches)  # horizontal, vertical
    height = 4.25*inches
    width = aspect*height
    syst['viewport'] = [cntr[0]-width/2.0,cntr[0]+width/2.0,cntr[1]-height/2.0,cntr[1]+height/2.0]
    fid.write(write_style.style2string(syst,landscape=1))
    fid.close()
    if n is None:
        winnum = gist.window(style=style,width=int(width*1.25/inches*_dpi),height=int(height*1.4/inches*_dpi))
        if winnum < 0:
            gist.window(style=style,width=int(width*1.25/inches*_dpi),height=int(height*1.4/inches*_dpi)) 
    else:
        gist.window(n,style=style,width=int(width*1.25/inches*_dpi),height=int(height*1.4/inches*_dpi))
        _current_style = style
    return

def full_page(win):
    gist.window(win,style=_current_style,width=int(_dpi*8.5),height=_dpi*11)

def _add_color(system, color, frame=0):
    try:
        system['ticks']['horiz']['tickStyle'] = {'color':color}
        system['ticks']['horiz']['gridStyle'] = {'color':color}
    except KeyError:
        system['ticks']['horiz'] = {}
        system['ticks']['horiz']['tickStyle'] = {'color':color}
        system['ticks']['horiz']['gridStyle'] = {'color':color}
    try:
        text = system['ticks']['horiz']['textStyle']
    except KeyError:
        system['ticks']['horiz']['textStyle'] = {}
    text = system['ticks']['horiz']['textStyle']
    text['color'] = color

    try:
        system['ticks']['vert']['tickStyle'] = {'color':color}
        system['ticks']['vert']['gridStyle'] = {'color':color}
    except KeyError:
        system['ticks']['vert'] = {}
        system['ticks']['vert']['tickStyle'] = {'color':color}
        system['ticks']['vert']['gridStyle'] = {'color':color}
    try:
        text = system['ticks']['vert']['textStyle']
    except KeyError:
        system['ticks']['vert']['textStyle'] = {}
    text = system['ticks']['vert']['textStyle']
    text['color'] = color

    system['ticks']['frame'] = frame
    system['ticks']['frameStyle'] = {'color':color}


def _chng_font(system, font, height):
    if height is None:
        height=14
    if font is None:
        font = 'helvetica'        
    num = write_style.tfont[font]
    system['ticks'] = { 
        'horiz':{
        'textStyle':{'font':num,
                     'height':height*points}
        },
        'vert':{
        'textStyle':{'font':num,
                     'height':height*points}
        }
    }
    return

def _remove_ticks(system):
    system['ticks'] = {'horiz': {'flags':0},
                       'vert': {'flags':0},
                       }
    return

plotframe = gist.plsys
import os
def subplot(Numy,Numx,win=0,pw=None,ph=None,hsep=100,vsep=100,color='black',frame=0,fontsize=8,font=None,ticks=1):
    # Use gist.plsys to change coordinate systems

    # all inputs (except fontsize) given as pixels, gist wants
    #  things in normalized device
    #  coordinate.  Window is brought up with center of window at
    #  center of 8.5 x 11 inch page: in landscape mode (5.25, 4.25)
    #  or at position (4.25,6.75) for portrait mode
    msg = 1
    if pw is None:
        pw = Numx*300
        msg = 0
    if ph is None:
        ph = Numy*300
        msg = 0
    maxwidth=int(os.environ.get('XPLT_MAXWIDTH'))
    maxheight=int(os.environ.get('XPLT_MAXHEIGHT'))


    printit = 0
    if ph > maxheight:
        ph = maxheight
        printit = 1
    if pw > maxwidth:
        pw = maxwidth
        printit = 1
        
    if _dpi != 100:
        fontsize = 12
    conv = inches *1.0 / _dpi  # multiply by this factor to convert pixels to
                              # NDC

    # Use landscape mode unless requested height is large
    land = 1
    maxw = 11*_dpi
    maxh = 8.5*_dpi
    if ph > (8.5*_dpi) and pw < (8.5*_dpi):
        land = 0
        maxh = 11*_dpi
        maxw = 8.5*_dpi
    
    if ph > maxh:
        ph = maxh
        printit=1
    if pw > maxw:
        pw = maxw
        printit=1

    if printit and msg:
        message = "Warning: Requested height and width too large.\n"
        message +="Changing to %d x %d" % (pw,ph)
        print message

    # Now we've got a suitable height and width

    if land:
        cntr = array([5.5,4.25])*_dpi  # landscape
    else:
        cntr = array([4.25,6.75])*_dpi  # portrait

    Yspace = ph/float(Numy)*conv
    Xspace = pw/float(Numx)*conv 

    hsep = hsep * conv
    vsep = vsep * conv
    ytop = (cntr[1]+ph/2.0)*conv
    xleft = (cntr[0]-pw/2.0)*conv

    if type(color) is types.StringType:
        color = _colornum[color]
    systems=[]
    ind = -1
    for nY in range(Numy):
        ystart = ytop - (nY+1)*Yspace
        for nX in range(Numx):
            xstart = xleft + nX*Xspace
            systems.append({})
            systems[-1]['viewport'] = [xstart+hsep/2.0,xstart+Xspace-hsep/2.0,ystart+vsep/2.0,ystart+Yspace-vsep/2.0]
            if font is not None or fontsize is not None:
                _chng_font(systems[-1],font,fontsize)
            if color != -3 or frame != 0:
                _add_color(systems[-1],color,frame=frame)
            if ticks != 1:
                _remove_ticks(systems[-1])
    _current_style=os.path.join(_user_path,"subplot%s.gs" % win)
    fid = open(_current_style,'w')
    fid.write(write_style.style2string(systems,landscape=land))
    fid.close()
    gist.winkill(win)
    gist.window(win,style=_current_style,width=int(pw),height=int(ph))

_dwidth=6*inches
_dheight=6*inches

import colorbar
def imagesc_cb(z,cmin=None,cmax=None,xryr=None,_style='default',
               zlabel=None,font='helvetica',fontsize=16,color='black',
               palette=None):
    """Plot an image on axes with a colorbar on the side.

    z -- The data
    cmin -- Value to map to lowest color in palette (min(z) if None)
    cmax -- Value to map to highest color in palette (max(z) if None)
    xryr -- (xmin, ymin, xmax, ymax) coordinates to print
            (0, 0, z.shape[1], z.shape[0]) if None
    _style -- A 'style-sheet' to use if desired (a default one will be used
              if 'default').  If None, then no style will be imposed.
    palette -- A string for a palette previously saved in a file (see write_palette)
               or an array specifying the red-green-blue values (2-d array N x 3) or
               gray-scale values (2-d array N x 1 or 1-d array).
    zlabel -- The label to attach to the colorbar (font, fontsize, and color
              match this).
    color -- The color to use for the ticks and frame.
    """    
    if xryr is None:
        xryr = (0,0,z.shape[1],z.shape[0])
        
    if not _hold:
        gist.fma()
    gist.animate(0)
    if _style is not None:
        if _style == 'default':
            _style=os.path.join(_user_path,"colorbar.gs")
            system = write_style.getsys(hticpos='below',vticpos='left',frame=1,color=color)
            fid = open(_style,'w')
            fid.write(write_style.style2string(system))
            fid.close()
        gist.window(style=_style)
        _current_style=_style
    if cmax is None:
        cmax = max(ravel(z))
    if cmin is None:
        cmin = min(ravel(z))        
    cmax = float(cmax)
    cmin = float(cmin)

    change_palette(palette)

    byteimage = gist.bytscl(z,cmin=cmin,cmax=cmax)
    gist.pli(byteimage,xryr[0],xryr[1],xryr[2],xryr[3])
    colorbar.color_bar(cmin,cmax,ncol=240,zlabel=zlabel,font=font,fontsize=fontsize,color=color)

def xlabel(text,color=None,font='helvetica',fontsize=16,deltax=0.0,deltay=0.0):
    """To get symbol font for the next character precede by !.  To get
    superscript enclose with ^^
    To get subscript enclose with _<text>_
    """
    global _textcolor
    if color is None:
        color = _textcolor
    else: 
        _textcolor = color
    if color is None:
        color = 'black'

    vp = gist.viewport()
    xmidpt = (vp[0] + vp[1])/2.0 + deltax
    y0 = vp[2] - 0.035 + deltay
    if text != "":
        gist.plt(text, xmidpt, y0, color=color,
                 font=font, justify="CT", height=fontsize)
    return xmidpt, y0


def ylabel(text,color=None,font='helvetica',fontsize=16,deltax=0.0,deltay=0.0):
    """To get symbol font for the next character precede by !.  To get
    superscript enclose with ^^
    To get subscript enclose with _<text>_
    """
    global _textcolor
    if color is None:
        color = _textcolor
    else:
        _textcolor = color
    if color is None:
        color = 'black'
    vp = gist.viewport()
    ymidpt = (vp[2] + vp[3])/2.0 + deltay
    x0 = vp[0] - 0.055 + deltax
    if text != "":
        gist.plt(text, x0, ymidpt, color=color,
                 font=font, justify="CB", height=fontsize, orient=1)
    return x0, ymidpt


def title(text,color=None,font='helvetica',fontsize=18,deltax=0.0,deltay=0.0):
    """Set title for plot.

    To get symbol font for the next character precede by !.  To get
    superscript enclose with ^^
    To get subscript enclose with _<text>_
    """
    global _textcolor
    if color is None:
        color = _textcolor
    else:
        _textcolor = color
    if color is None:
        color = 'black'
        
    vp = gist.viewport()
    xmidpt = (vp[0] + vp[1])/2.0 + deltax
    if text != "":
        gist.plt(text,xmidpt,vp[3] + 0.02 + deltay, font=font, justify='CB',
                 height=fontsize, color=color)

def title3(text,color=None,font='helvetica',fontsize=18,deltax=0.0,deltay=0.0):

    global _textcolor
    if color is None:
        color = _textcolor
    else:
        _textcolor = color
    if color is None:
        color = 'black'

    vp = gist.viewport()
    xmidpt = (vp[0] + vp[1])/2.0 + deltax
    if text != "":
        gist.plt(text,xmidpt,vp[3]-0.05-deltay, font=font, justify='CB',
                 height=fontsize, color=color)

def stem(m, y, linetype='b-', mtype='mo', shift=0.013):
    y0 = Numeric.zeros(len(y),y.typecode())
    y1 = y
    x0 = m
    x1 = m
    try:
        override = 1
        savesys = gist.plsys(2)
        gist.plsys(savesys)
    except:
        override = 0
    if not (_hold or override):
        gist.fma()
    thetype,thecolor,themarker,tomark = _parse_type_arg(linetype,0)
    lcolor = thecolor
    gist.pldj(x0, y0, x1, y1, color=thecolor, type=thetype)
    thetype,thecolor,themarker,tomark = _parse_type_arg(mtype,0)
    if themarker not in ['o','x','.','*']:
        themarker = 'o'
    y = where(scipy.isfinite(y),y,0)
    gist.plg(y,m,color=thecolor,marker=themarker,type='none')
    gist.plg(Numeric.zeros(len(m)),m,color=lcolor,marks=0)
    gist.limits()
    lims = gist.limits()
    newlims = [None]*4
    vp = gist.viewport()
    factor1 = vp[1] - vp[0]
    factor2 = vp[3] - vp[2]
    cfactx = factor1 / (lims[1] - lims[0])
    cfacty = factor2 / (lims[3] - lims[2])
    d1 = shift / cfactx
    d2 = shift / cfacty
    newlims[0] = lims[0] - d1
    newlims[1] = lims[1] + d1
    newlims[2] = lims[2] - d2
    newlims[3] = lims[3] + d2
    gist.limits(*newlims)
    return


def makeleg(leg,pos,lenx,dd,theight=12):
    # Place legend
    x0,y0 = pos
    dx,dy = dd
    for k in range(len(leg['txt'])):
        gist.plg([y0+k*dy]*2,[x0,x0+lenx],type=leg['sym'][k][1],marks=0)
        if leg['sym'][k][0] is not None:
            gist.plg([y0+k*dy]*2,[x0,x0+lenx],type='none',marks=1,marker=leg['sym'][k][0])
        if leg['txt'][k] != "":
            gist.plt(leg['txt'][k],x0+lenx+dx,y0+k*dy,height=theight,tosys=1,justify='LH') 
    return

def twoplane(DATA,slice1,slice2,dx=[1,1,1],cmin=None,cmax=None,xb=None,xe=None,
             xlab="",ylab="",zlab="",clab="",titl="",
             totalheight=0.5,space=0.02, medfilt=5,
             font='helvetica',fontsize=16,color='black',lcolor='white',
             fcolor='black',  cb=1, line=1, palette=None):
    """ Visualize a 3d volume as a two connected slices.  The slices are
    given in the 2-tuple slice1 and slice2.

    These give the dimension and corresponding slice numbers to plot.
    The unchosen slice is the common dimension in the images.

    twoplane(img3d,(0,12),(2,60)) plots two images with a common "x"-axis
    as the first dimension.  The lower plot is img3d[12,:,:] with a line
    through row 60 corresponding to the slice transpose(img3d[:,:,60])
    plotted above this first plot.
    """
    if xb is None:
        xb = [0,0,0]
    if xe is None:
        xe = DATA.shape
    # get two image slices 
    # make special style file so that pixels are square
    getdx = array([1,1,1])
    imgsl1 = [slice(None,None),slice(None,None),slice(None,None)]
    imgsl1[slice1[0]] = slice1[1]
    img1 = DATA[imgsl1]
    getdx1 = getdx.__copy__()
    getdx1[slice1[0]] = 0
    dx1 = compress(getdx1,dx)
    xb1 = compress(getdx1,xb)
    xe1 = compress(getdx1,xe)

    imgsl2 = [slice(None,None),slice(None,None),slice(None,None)]
    imgsl2[slice2[0]] = slice2[1]
    img2 = DATA[imgsl2]
    getdx2 = getdx.__copy__()
    getdx2[slice2[0]] = 0
    dx2 = compress(getdx2,dx)
    xb2 = compress(getdx2,xb)
    xe2 = compress(getdx2,xe)


    if (slice1[0] == slice2[0]):
        raise ValueError, "Same slice dimension.."

    for k in range(3):
        if k not in [slice1[0],slice2[0]]:
            samedim = k
            break
    if samedim == 2:
        pass
    elif samedim == 1:
        if samedim > slice1[0]:
            img1 = transpose(img1)
            dx1 = dx1[::-1]
            xb1 = xb1[::-1]
            xe1 = xe1[::-1]
        if samedim > slice2[0]:
            img2 = transpose(img2)
            dx2 = dx2[::-1]
            xb2 = xb2[::-1]
            xe2 = xe2[::-1]
    else:
        img1 = transpose(img1)
        dx1 = dx1[::-1]
        xb1 = xb1[::-1]
        xe1 = xe1[::-1]
        img2 = transpose(img2)
        dx2 = dx2[::-1]
        xb2 = xb2[::-1]
        xe2 = xe2[::-1]
        
        

    assert(img1.shape[1] == img2.shape[1])
    units = totalheight - space
    totaldist = img1.shape[0]*dx1[0] + img2.shape[0]*dx2[0]
    convfactor = units / float(totaldist)
    height1 = img1.shape[0]*dx1[0] * convfactor
    xwidth = img1.shape[1]*dx1[1]*convfactor
    if xwidth > 0.6:
        rescale = 0.6 / xwidth
        xwidth = rescale * xwidth
        height1 = rescale * height1
        totalheight = totalheight * rescale
        print xwidth, height1
    else:
        print xwidth
    ystart = 0.5 - totalheight / 2
    ypos1 = [ystart, ystart+height1]
    ypos2 = [ystart+height1+space,ystart+totalheight]    
    xpos = [0.395-xwidth/2.0, 0.395+xwidth/2.0]

    systems = []
    system = write_style.getsys(hticpos='', vticpos='left')
    system['viewport'] = [xpos[0],xpos[1],ypos2[0],ypos2[1]]
    if fcolor not in ['black',None]:
        _add_color(system, _colornum[color])
    systems.append(system)
    system = write_style.getsys(hticpos='below', vticpos='left')
    system['viewport'] = [xpos[0],xpos[1],ypos1[0],ypos1[1]]
    if fcolor not in ['black',None]:
        _add_color(system, _colornum[color])
    systems.append(system)

    the_style = os.path.join(_user_path,"two-plane.gs")
    write_style.writestyle(the_style,systems)

    gist.window(style=the_style)
    _current_style = the_style
    change_palette(palette)
    gist.plsys(1)
    if medfilt > 1:
        img1 = signal.medfilt(img1,[medfilt,medfilt])
        img2 = signal.medfilt(img2,[medfilt,medfilt])
    if cmax is None:
        cmax = max(max(ravel(img1)),max(ravel(img2)))
    if cmin is None:
        cmin = min(min(ravel(img1)),min(ravel(img2)))
    cmax = float(cmax)
    cmin = float(cmin)
    byteimage = gist.bytscl(img2,cmin=cmin,cmax=cmax)
    gist.pli(byteimage,xb2[1],xb2[0],xe2[1],xe2[0])
    ylabel(zlab,color=color)
    if titl != "":
        title(titl,color=color)
    if line:
        xstart = xb2[1]
        xstop = xe2[1]
        yval = slice1[1]*(xe2[0] - xb2[0])/(img2.shape[0]) + xb2[0]
        gist.pldj([xstart],[yval],[xstop],[yval],type='dash',width=2,color='white')


    gist.plsys(2)
    ylabel(ylab,color=color)
    xlabel(xlab,color=color)
    byteimage = gist.bytscl(img1,cmin=cmin,cmax=cmax)
    gist.pli(byteimage,xb1[1],xb1[0],xe1[1],xe1[0])
    if line:
        xstart = xb1[1]
        xstop = xe1[1]
        yval = slice2[1]*(xe1[0] - xb1[0])/(img1.shape[0]) + xb1[0]
        gist.pldj([xstart],[yval],[xstop],[yval],type='dash',width=2,color='white')

    if cb:
        colorbar.color_bar(cmin,cmax,ncol=240,zlabel=clab,font=font,fontsize=fontsize,color=color,ymin=ystart,ymax=ystart+totalheight,xmin0=xpos[1]+0.02,xmax0=xpos[1]+0.04) 

def surf(z,x=None,y=None,win=None,shade=0,edges=1,edge_color="black",phi=-45.0,
         theta=30.0,zscale=1.0,palette=None,gnomon=0):
    """Plot a three-dimensional wire-frame (surface): z=f(x,y)
    """
    if win is None:
        pl3d.window3()
    else:
        pl3d.window3(win)
    pl3d.set_draw3_(0)
    phi0 = phi*pi/180.0
    theta0 = theta*pi/180.0
    pl3d.orient3(phi=phi0,theta=theta0)
    pl3d.light3()
    change_palette(palette)
    sz = scipy_base.shape(z)
    if len(sz) != 2:
        raise ValueError, "Input must be a 2-d array --- a surface."
    N,M = sz
    if x is None:
        x = arange(0,N)
    if y is None:
        y = arange(0,M)
    x = scipy_base.squeeze(x)
    y = scipy_base.squeeze(y)
    if (len(scipy_base.shape(x)) == 1):
        x = x[:,NewAxis]*ones((1,M))
    if (len(scipy_base.shape(y)) == 1):
        y = ones((N,1))*y[NewAxis,:]
    plwf.plwf(z,y,x,shade=shade,edges=edges,ecolor=edge_color,scale=zscale)
    lims = pl3d.draw3(1)
    gist.limits(lims[0],lims[1],lims[2],lims[3])
    pl3d.gnomon(gnomon)

def expand_limits(xpcnt,ypcnt=None):
    """Expand the limits by a certain percentage.
    """
    if ypcnt is None:
        ypcnt = xpcnt
    if xpcnt > 1:
        xpcnt = xpcnt / 100.0
    if ypcnt > 1:
        ypcnt = ypcnt / 100.0
    xmin, xmax, ymin, ymax, flag = gist.limits()
    dx = (xmax-xmin)*xpcnt/2.0
    dy = (ymax-ymin)*ypcnt/2.0
    gist.limits(xmin-dx,xmax+dx,ymin-dy,ymax+dy)

def axes(type='b|'):
    vals = gist.limits()
    v0,v1,v2,v3 = vals[:4]
    x0 = scipy_base.r_[v0:v1:5j]
    y0 = 5*[0]
    x1 = 5*[0]
    y1 = scipy_base.r_[v2:v3:5j]
    plot(x0,y0,type,x1,y1,type,hold=1)
    gist.limits(v0,v1,v2,v3)
    

def bode(w,H,win=0,frame=0,lcolor='blue',color='black',tcolor='black',freq='rad'):
    """Plot a bode plot of the transfer function H as a function of w.
    """
    if freq == 'Hz':
        w = w /2.0 / pi
    subplot(2,1,win,hsep=120,frame=frame,color=color)
    gist.plsys(1)
    gist.plg(20*scipy.log10(abs(H)),w,type='solid',color=lcolor,marks=0)
    gist.logxy(1,0)
    gist.gridxy(1,1)
    if freq == 'Hz':
        xlabel('Frequency (Hz)',color=tcolor,deltay=-0.005)
    else:
        xlabel('Frequency (rad/s)',color=tcolor,deltay=-0.005)         
    ylabel('Magnitude (dB)',color=tcolor,deltax=-0.005)
    title("Bode Plot",color=tcolor)
    gist.plsys(2)
    gist.plg(180/pi*scipy.unwrap(MLab.angle(H)),w,type='solid',color=lcolor,marks=0)
    gist.logxy(1,0)
    gist.gridxy(1,1)
    if freq == 'Hz':
        xlabel('Frequency (Hz)',color=tcolor,deltay=-0.005)
    else:
        xlabel('Frequency (rad/s)',color=tcolor,deltay=-0.005)         
    ylabel('Phase (deg.)',color=tcolor,deltax=-0.005)
    

def addtext(txt,xy=None,fontsize=16,font='helvetica',color='black',
            orient=0,justify='LA',tosys=0):
    if xy is None:
        result = gist.mouse(0,0,"Click on point for lower left starting position")
        if result is None:
            raise ValueError, "Invalid point entered."        
        x,y = result[4],result[5]
        tosys = 0
    else:
        x,y = xy
    if isinstance(txt, type('')) and txt != '':
        gist.plt(txt, x, y, tosys=tosys, justify=justify,
                 height=fontsize,font=font, color=color, orient=orient)
    return

    
