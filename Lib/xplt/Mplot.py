# Author: Travis Oliphant

import gist
import pl3d, plwf
import Numeric
from Numeric import ravel, reshape, repeat, arange, transpose, compress
import MLab
from MLab import pi, cos, sin, arctan2, array
import types
import write_style
from write_style import inches, points
import scipy
from scipy.signal import medfilt

_hold = 0

try:
    import Scientific.Statistics.Histogram
    SSH = Scientific.Statistics.Histogram
    def histogram(data,numbins=80,range=None):
        """Plot a histogram.
        """

        h = SSH.Histogram(data,numbins,range)
        plot(h[:,0],h[:,1])
        return h
except ImportError:
    try:
        import Statistics
        SSH = Statistics
        def histogram(data,numbins=80,range=None):
            """Plot a histogram.
            """
            h = SSH.histogram(data,numbins,range)
            plot(h[:,0],h[:,1])
            return h        
    except ImportError:
        histogram = scipy.histogram

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
    if not _hold:
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
    if _hold:
        pass
    else:
        gist.fma()
    gist.plg(y,x,color=_colors[ptcolor],marker=_markers[pttype],type='none')
    gist.pldj(x,yb,x,ye,color=_colors[linecolor],type=_types[linetype])
    viewp = gist.viewport()
    plotlims = gist.limits()
    conv_factorx = (viewp[1] - viewp[0]) / (plotlims[1]-plotlims[0])
    conv_factory = (viewp[3] - viewp[2]) / (plotlims[3]-plotlims[2])    
    width = fac*(x[2]-x[1])
    x0 = x-width/2.0
    x1 = x+width/2.0
    gist.pldj(x0,ye,x1,ye,color=_colors[linecolor],type=_types[linetype])
    gist.pldj(x0,yb,x1,yb,color=_colors[linecolor],type=_types[linetype])
    return

def legend(text,linetypes=None,lleft=None,color='black',tfont='helvetica',fontsize=14,nobox=0):
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
    viewp = gist.viewport()
    width = (viewp[1] - viewp[0]) / 10.0;
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
    dy = fontsize*points*1.15
    deltay = fontsize*points / 2.8
    deltax = fontsize*points / 2.8
    ypos = lly + deltay;
    if linetypes is None:
        linetypes = _GLOBAL_LINE_TYPES[:]  # copy them out
    gist.plsys(0)
    savehold = _hold
    _hold = 1
    for k in range(len(text)):
        plot(legarr,ypos*legy,linetypes[k])
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
    if "hold" in keywds.keys():
        _hold = keywds['hold']
    if _hold or override:
        pass
    else:
        gist.fma()
    gist.animate(0)
    nargs = len(args)
    if nargs == 0:
        y = x
        x = Numeric.arange(0,len(y))
        gist.plg(y,x,type='solid',color='blue',marks=0)
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
        gist.plg(y,x,type=thetype,color=thecolor,marker=themarker,marks=tomark)

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
    assert(len(y.shape)==2)
    assert(len(x)==y.shape[axis])
    otheraxis = (1+axis) % 2
    sliceobj = [slice(None)]*2
    if not _hold:
        gist.fma()
    for k in range(y.shape[otheraxis]):
        thiscolor = _colors[_corder[k % len(_corder)]] 
        sliceobj[otheraxis] = k
        gist.plg(y[sliceobj],x,type='solid',color=thiscolor,marks=0)


def addbox(x0,y0,x1,y1,color='black',width=1,type='-'):
    if not isinstance(color,types.IntType):
        color = _colornum[color]
    wordtype = _types[type]
    gist.pldj([x0,x1,x1,x0],[y0,y0,y1,y1],[x1,x1,x0,x0],[y0,y1,y1,y0],color=color,type=wordtype,width=width)

def imagesc(z,cmin=None,cmax=None,xryr=None,_style=None,mystyle=0,
            palette=None):
    if xryr is None:
        xryr = (0,0,z.shape[1],z.shape[0])
    if not _hold:
        gist.fma()
    gist.animate(0)
    if _style is None and mystyle==0:
        _style='/tmp/image.gs'
        system = write_style.getsys(hticpos='left',frame=1)
        fid = open(_style,'w')
        fid.write(write_style.style2string(system))
        fid.close()
    gist.window(style=_style)
    if cmax is None:
        cmax = max(ravel(z))
    if cmin is None:
        cmin = min(ravel(z))
    cmax = float(cmax)
    cmin = float(cmin)
    byteimage = gist.bytscl(z,cmin=cmin,cmax=cmax)
    if palette is not None:
        try:
            gist.palette('%d.gp' % palette)
        except IOError:
            raise ValueError, "Palette %d not found."
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

def figure(n=None, style='/tmp/currstyle.gs', color=-2, frame=1, labelsize=14, labelfont='helvetica'):
    fid = open(style,'w')
    fid.write(write_style.style2string(write_style.getsys(color=color,frame=frame,labelsize=labelsize,font=labelfont)))
    fid.close()
    if n is None:
        gist.window(style=style)
    else:
        gist.window(n,style=style)
    return

def full_page(win,dpi=75):
    gist.window(win,style=_current_style,width=int(dpi*8.5),height=dpi*11,dpi=dpi)

plotframe = gist.plsys
def subplot(Numy,Numx,win=0,lm=0*inches,rm=0*inches,tm=0*inches,bm=0*inches,ph=11*inches,pw=8.5*inches,dpi=75,ls=0.75*inches,rs=0.75*inches,ts=0.75*inches,bs=0.75*inches):
    # Use gist.plsys to change coordinate systems 
    systems=[]
    ind = -1
    Yspace = (ph-bm-tm)/float(Numy)
    Xspace = (pw-rm-lm)/float(Numx)
    for nY in range(Numy):
        ystart = (ph-tm) - (nY+1)*Yspace + bs
        for nX in range(Numx):
            xstart = lm + nX*Xspace + ls
            systems.append({})
            systems[-1]['viewport'] = [xstart,xstart+Xspace-(ls+rs),ystart,ystart+Yspace-(ts+bs)]
    _current_style='/tmp/subplot%s.gs' % win
    fid = open(_current_style,'w')
    fid.write(write_style.style2string(systems))
    fid.close()
    gist.winkill(win)
    gist.window(win,style=_current_style,width=int(dpi*8.5),height=dpi*11,dpi=dpi)

_dwidth=6*inches
_dheight=6*inches

import colorbar
def imagesc_cb(z,cmin=None,cmax=None,xryr=None,_style=None,mystyle=0,
               zlabel=None,font='helvetica',fontsize=16,color='black'):
    if xryr is None:
        xryr = (0,0,z.shape[1],z.shape[0])
    if not _hold:
        gist.fma()
    gist.animate(0)
    if _style is None and mystyle==0:
        _style='/tmp/colorbar.gs'
        system = write_style.getsys(hticpos='left',frame=1)
        fid = open(_style,'w')
        fid.write(write_style.style2string(system))
        fid.close()
    gist.window(style=_style)
    if cmax is None:
        cmax = max(ravel(z))
    if cmin is None:
        cmin = min(ravel(z))        
    cmax = float(cmax)
    cmin = float(cmin)

    byteimage = gist.bytscl(z,cmin=cmin,cmax=cmax)
    gist.pli(byteimage,xryr[0],xryr[1],xryr[2],xryr[3])
    colorbar.color_bar(cmin,cmax,ncol=240,zlabel=zlabel,font=font,fontsize=fontsize,color=color)

def xlabel(text,color='black',font='helvetica',fontsize=16,deltax=0.0,deltay=0.0):
    vp = gist.viewport()
    xmidpt = (vp[0] + vp[1])/2.0 + deltax
    y0 = vp[2] - 0.040 + deltay
    if text != "":
        gist.plt(text, xmidpt, y0, color=color,
                 font=font, justify="CT", height=fontsize)
    return xmidpt, y0


def ylabel(text,color='black',font='helvetica',fontsize=16,deltax=0.0,deltay=0.0):
    vp = gist.viewport()
    ymidpt = (vp[2] + vp[3])/2.0 + deltay
    x0 = vp[0] - 0.040 + deltax
    if text != "":
        gist.plt(text, x0, ymidpt, color=color,
                 font=font, justify="CB", height=fontsize, orient=1)
    return x0, ymidpt


def title(text,color='black',font='helvetica',fontsize=18,deltax=0.0,deltay=0.0):
    vp = gist.viewport()
    xmidpt = (vp[0] + vp[1])/2.0 + deltax
    if text != "":
        gist.plt(text,xmidpt,vp[3] + 0.02 + deltay, font=font, justify='CB',
                 height=fontsize, color=color)

def title3(text,color='black',font='helvetica',fontsize=18,deltax=0.0,deltay=0.0):
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
             cb=1, line=1):
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
    systems.append(system)
    system = write_style.getsys(hticpos='below', vticpos='left')
    system['viewport'] = [xpos[0],xpos[1],ypos1[0],ypos1[1]]
    systems.append(system)

    write_style.writestyle("/tmp/two-plane.gs",systems)

    gist.window(style='/tmp/two-plane.gs')
    gist.plsys(1)
    if medfilt > 1:
        img1 = medfiltND(img1,[medfilt,medfilt])
        img2 = medfiltND(img2,[medfilt,medfilt])
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
def plot3(x,y,z,win=None,shade=0,edges=1,edge_color="black",phi=-45,theta=30,
          zscale=1.0,palette=None,gnomon=0):
    """Plot a three-dimensional wire-frame (surface): z=f(x,y)
    """
    if win is None:
        pl3d.window3()
    else:
        pl3d.window3(win)
    pl3d.set_draw3_(0)
    pl3d.orient3(phi=phi*pi/180,theta=theta*pi/180)
    pl3d.light3()
    if palette is not None:
        try:
            gist.palette('%s.gp' % palette)
        except:
            raise ValueError, "%s palette not available." % palette
    plwf.plwf(z,y,x,shade=shade,edges=edges,ecolor=edge_color,scale=zscale)
    [xmin,xmax,ymin,ymax] = pl3d.draw3(1)
    gist.limits(xmin,xmax,ymin,ymax)
    pl3d.gnomon(gnomon)
    


















