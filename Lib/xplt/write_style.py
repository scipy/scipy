import types
points = 0.0013000
inches = 72.27*points
def_linewidth = 0.5*points
def_msize = 10.0*points
cm = inches/2.54

ltype = {'none':0,
         'solid':1,
         'dash':2,
         'dot':3,
         'dashdot':4,
         'dashdotdot':5}
mtype = {'point':1,
         'plus':2,
         'asterisk':3,
         'circle':4,
         'cross':5}
tfont = {'bold':1,
         'italic':2,
         'courier':0,
         'times':4,
         'helvetica':8,
         'symbol':12,
         'newcentury':16}
fontnames = ['courier', 'times', 'helvetica', 'symbol', 'newcentury']
fontstyles = ['bold', 'italic']
tdirec = {'right':0,
          'up':1,
          'left':2,
          'down':3}
thalign = {'normal':0,
           'left':1,
           'center':2,
           'right':3}
tvalign = {'normal':0,
           'top':1,
           'cap':2,
           'half':3,
           'base':4,
           'bottom':5}
_ticks = {'none':0,
         'left':0x001,
         'below':0x001,
         'right':0x002,
         'above':0x002,
         'center':0x004,
         'in':0x008,
         'out':0x010,
         'leftlabel':0x020,
         'belowlabel':0x020,
         'rightlabel':0x040,
         'abovelabel':0x040,
         'fullgrid':0x080,
         'origin':0x100}
colors = {'bg':-1,
          'fg':-2,
          'black':-3,
          'white':-4,
          'red':-5,
          'green':-6,
          'blue':-7,
          'cyan':-8,
          'magenta':-9,
          'yellow':-10}

defsys = {'legend':0,
          'viewport':[0.19,0.60,0.44,0.85],
          'ticks':{
              'horiz':{
                  'nMajor':7.5,
                  'nMinor':50.0,
                  'logAdjMajor':1.2,
                  'logAdjMinor':1.2,
                  'nDigits':3,
                  'gridLevel':1,
                  'flags':0x033,
                  'tickOff':0.0007,
                  'labelOff':0.0182,
                  'tickLen':[0.0143,0.0091,0.0052,0.0026,0.0013],
                  'tickStyle':{
                      'color':-2,
                      'type':1,
                      'width':1.0},
                  'gridStyle':{
                      'color':-2,
                      'type':3,
                      'width':1.0},
                  'textStyle':{
                      'color':-2,
                      'font':0x08,
                      'height':0.0182,
                      'orient':0,
                      'alignH':0,
                      'alignV':0,
                      'opaque':0},
                  'xOver':0.395,
                  'yOver':0.370
                  },
              'vert':{
                  'nMajor':7.5,
                  'nMinor':50.0,
                  'logAdjMajor':1.2,
                  'logAdjMinor':1.2,
                  'nDigits':4,
                  'gridLevel':1,
                  'flags':0x033,
                  'tickOff':0.0007,
                  'labelOff':0.0182,
                  'tickLen':[0.0143,0.0091,0.0052,0.0026,0.0013],
                  'tickStyle':{
                      'color':-2,
                      'type':1,
                      'width':1.0},
                  'gridStyle':{
                      'color':-2,
                      'type':3,
                      'width':1.0},
                  'textStyle':{
                      'color':-2,
                      'font':0x08,
                      'height':0.0182,
                      'orient':0,
                      'alignH':0,
                      'alignV':0,
                      'opaque':0},
                  'xOver':0.150,
                  'yOver':0.370
                  },
              'frame':0,
              'frameStyle':{
                  'color':-2,
                  'type':1,
                  'width':1.0}
              }
          }

def getsys(color=-2,frame=1,labelsize=14, font='helvetica', ticks='solid', hticpos='left right', vticpos='below above'):
    if type(color) != types.IntType:
        color = colors[color]
    labels = labelsize * points
    fontlist = string.split(string.lower(font))
    found = 0
    for k in range(len(fontlist)):
        if fontlist[k] in fontnames:
            found = 1
            break

    if not found:
        print "Unknown font: %s" % font
        return

    Nfont = tfont[fontlist[k]]
    for k in range(len(fontstyles)):
        if fontstyles[k] in fontlist:
            Nfont = Nfont + tfont[fontstyles[k]]

    ticksty = ltype[ticks]
    if ticksty == 0:
        hflags = 0
        vflags = 0
    else:
        hflags = vflags = 0x030
        for pos in string.split(hticpos):
            try:
                hflags = hflags + _ticks[pos]
            except KeyError:
                pass
        for pos in string.split(vticpos):
            try:
                vflags = vflags + _ticks[pos]
            except KeyError:
                pass

    newsys = {'legend':0,
              'viewport':[0.19,0.60,0.44,0.85],
              'ticks':{
                  'horiz':{
                      'nMajor':7.5,
                      'nMinor':50.0,
                      'logAdjMajor':1.2,
                      'logAdjMinor':1.2,
                      'nDigits':3,
                      'gridLevel':1,
                      'flags':hflags,
                      'tickOff':0.0007,
                      'labelOff':0.0182,
                      'tickLen':[0.0143,0.0091,0.0052,0.0026,0.0013],
                      'tickStyle':{
                          'color':color,
                          'type':ticksty,
                          'width':1.0},
                      'gridStyle':{
                          'color':color,
                          'type':3,
                          'width':1.0},
                      'textStyle':{
                          'color':color,
                          'font':Nfont,
                          'height':labels,
                          'orient':0,
                          'alignH':0,
                          'alignV':0,
                          'opaque':0},
                      'xOver':0.395,
                      'yOver':0.370
                      },
                  'vert':{
                      'nMajor':7.5,
                      'nMinor':50.0,
                      'logAdjMajor':1.2,
                      'logAdjMinor':1.2,
                      'nDigits':4,
                      'gridLevel':1,
                      'flags':vflags,
                      'tickOff':0.0007,
                      'labelOff':0.0182,
                      'tickLen':[0.0143,0.0091,0.0052,0.0026,0.0013],
                      'tickStyle':{
                          'color':color,
                          'type':ticksty,
                          'width':1.0},
                      'gridStyle':{
                          'color':color,
                          'type':3,
                          'width':1.0},
                      'textStyle':{
                          'color':color,
                          'font':Nfont,
                          'height':labels,
                          'orient':0,
                          'alignH':0,
                          'alignV':0,
                          'opaque':0},
                      'xOver':0.150,
                      'yOver':0.370
                      },
                  'frame':frame,
                  'frameStyle':{
                      'color':color,
                      'type':1,
                      'width':1.0}
                  }
              }
    return newsys
    
import string, types

def sys2string(system,level=0):
    retstr=""
    spaces = 4*level
    for key in system.keys():
        keytype = type(system[key])
        if keytype is types.DictType:
            retstr = retstr+"%s%s={\n%s" % (' '*spaces,key,sys2string(system[key],level+1))
        elif keytype is types.ListType:
            retstr = retstr+"%s%s={%s},\n" % (' '*spaces,key,string.join(map(str,system[key]),','))
        else:
            retstr = retstr+"%s%s=%s,\n" % (' '*spaces,key,system[key])
    return retstr[:-2]+"\n%s},\n" % (' '*spaces)
            
def style2string(systemslist, landscape=0):
    retstr = "landscape = %d\n" % landscape
    num = 1
    if type(systemslist) is types.ListType:
        num = len(systemslist)
        if type(systemslist[0]) is not types.DictType:
            raise TypeError, "dict2string: first argument should be a dicitonary or a list"
    elif type(systemslist) is not types.DictType:
        raise TypeError, "dict2string: first argument should be a dictionary or a list"
    else:
        systemslist = [systemslist]

    if num > 1:
        retstr = "%s\ndefault = {\n%s}\n" % (retstr, sys2string(systemslist[0],level=1)[:-3])
        retstr = "%s\nsystem = { legend=0 }\n" % retstr
        for k in range(1,num):
            retstr = "%s\nsystem = {\n%s}\n" % (retstr, sys2string(systemslist[k],level=1)[:-3])
    else:
        retstr = "%s\nsystem = {\n%s}\n" % (retstr, sys2string(systemslist[0],level=1)[:-3])

    return retstr

def writestyle(name, systemslist, landscape=0):
    fid = open(name,'w')
    fid.write(style2string(systemslist,landscape))
    fid.close()
    return



