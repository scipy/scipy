# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.
from Numeric import *
from fastumath import *
import sys, os	# To be sure expand_path has posixpath and we have sys.path
try:
    from gistC import *
except ImportError:
    from scipy.xplt.gistC import *
from helpmod import help as ghelp
from shapetest import *
from arrayfns import *

# Parameters used by pltitle and xytitles
pltitle_height= 18;
pltitle_font= "helvetica";

# Parameters used by plmk and plmk_default
_plmk_count = 0
_plmk_msize = _plmk_color = _plmk_width = None

# delete the current graphics window, or graphics window N (0-7).
def winkill(*N):
  if N:
    window(N[0], display="", hcp="")
  else:
    window(display="", hcp="")

# Plot TITLE centered above the coordinate system for any of the
# standard Gist styles.  You will need to customize this for other
# plot styles.
def pltitle(title):
  vp = viewport()
  xmidpt = (vp[0] + vp[1])/2.0
  plt( title, xmidpt, vp[3] + 0.02,
       font=pltitle_font, justify="CB", height=pltitle_height)

# Set limits on the y-axis
def ylimits(ymin='u',ymax='u'): limits('u','u',ymin,ymax)

# Return the 1-origin zone index for the point clicked in
# for the default mesh, or for the mesh (X,Y) (region array IREG).
# Number of args can be 0, 2, or 3.
def moush(*arg):
  narg = len(arg)
  if narg == 3: # (y, x, ireg)
    xy = mouse (-1, 0, "<Click mouse in mesh>")
    if xy is None: return None
    return mesh_loc (xy[1], xy[0], arg[0], arg[1], arg[2]);
  elif narg == 2: # (y, x)
    xy = mouse (-1, 0, "<Click mouse in mesh>")
    if xy is None: return None
    return mesh_loc (xy[1], xy[0], arg[0], arg[1]);
  elif narg == 0: # ()
    xy = mouse (-1, 0, "<Click mouse in mesh>")
    if xy is None: return None
    return mesh_loc (xy[1], xy[0]);
  else:
    print "Moush takes 0, 2, or 3 args: ( [ y, x [ , ireg ] ] )"
    return None

# Create an encapsulated PostScript file.  Requires Ghostscript and its
# associated ps2epsi utility.
def eps(name, nops2epsi=0, pdf=0):
  import os
  totalname = name
  apath, basename = os.path.split(name)
  name = name + ".ps"
  window (hcp = name, dump = 1, legends = 0)
  hcp ()
  window (hcp="")
  if not nops2epsi:
      os.system ("ps2epsi " + name)
      os.system ("rm " + name)
      os.system ("mv %s.epsi %s.eps" % (basename, totalname))
  else:
      os.system ("mv %s.ps %s.eps" % (totalname, totalname))
  if pdf:
      os.system("epstopdf %s.eps" % totalname)

def xytitles(xtitle = "", ytitle = "", delta = (0.,0.)):
  vp = viewport()
  xmidpt = (vp[0] + vp[1])/2.0
  ymidpt = (vp[2] + vp[3])/2.0
  if len(xtitle) > 0:
    plt(xtitle, xmidpt, vp[2] - 0.050 + delta[1],
        font=pltitle_font, justify="CT", height=pltitle_height)
  if len(ytitle) > 0:
    plt(ytitle, vp[0] - 0.050 + delta[0], ymidpt,
        font=pltitle_font, justify="CB", height=pltitle_height, orient=1)

def _spanz(lb,ub,n):
  if n < 3: raise ValueError, '3rd arg must be at least 3'
  c = 0.5*(ub - lb)/(n - 1.0)
  b = lb + c
  a = (ub - c - b)/(n - 2.0)
  return map(lambda x,A=a,B=b: A*x + B, arange(n-1))

# predefined markers: square, +, delta, circle, diamond, x, grad
_seq = _spanz(-pi,pi,37)
_plmk_markers = (
  array([[-1,1,1,-1],[-1,-1,1,1]])*.007,
  array([[-4,-1,-1,1,1,4,4,1,1,-1,-1,-4],
    [-1,-1,-4,-4,-1,-1,1,1,4,4,1,1]])*.007/sqrt(7),
  array([[-sqrt(3),sqrt(3),0],[-1,-1,2]])*.007/sqrt(.75*sqrt(3)),
  array([cos(_seq),sin(_seq)])*.007/(pi/4.),
  array([[-1,0,1,0],[0,-1,0,1]])*.007*sqrt(2),
  array([[-1,-2.5,-1.5,0,1.5,2.5,1,2.5,1.5,0,-1.5,-2.5],
    [0,-1.5,-2.5,-1,-2.5,-1.5,0,1.5,2.5,1,2.5,1.5]])*.007*sqrt(2)/sqrt(7),
  array([[0,sqrt(3),-sqrt(3)],[-2,1,1]])*.007/sqrt(.75*sqrt(3))
  )
del(_seq)

def plmk(y,x=None,marker=None,width=None,color=None,msize=None):
  global _plmk_count
  global _plmk_color, _plmk_width, _plmk_msize
  color_dict = { 'bg':-1, 'fg':-2, 'black':-3, 'white':-4, 'red':-5,
    'green':-6, 'blue':-7, 'cyan':-8, 'magenta':-9, 'yellow':-10 }

  z = None

  if marker is None:
    marker = _plmk_markers[(_plmk_count)%7]
    _plmk_count = _plmk_count + 1
  elif type(marker) == type(0):
    marker = _plmk_markers[marker-1];

  xm = marker[0]
  ym = marker[1]
  if not msize: msize = _plmk_msize;
  if msize:
    xm = xm * msize;
    ym = ym * msize;

  if not color: color = _plmk_color;
  ecolor = color;
  if type(color) == type(""):
    color = color_dict[color];
  
  if not width: width = _plmk_width;
  if width >= 10:
    if not color:
      color = ecolor = -2 # magic number for "fg"
    z = ones(1+len(y)) * color
    z = z.astype(UnsignedInt8) # convert array to type <unsigned char>
    width = None

  n = ones(1+len(y));
  n[0] = len(ym);
  if not x: x = 1 + arange(len(y));
  plfp( z, concatenate((ym,y)), concatenate((xm,x)),n,
        edges=1, ewidth=width, ecolor=ecolor)

# Set default color, msize, and width values for plmk.  Use
# width=10 to get solid fills.  With no parameters, plmk_default
# restores the initial default values.
def plmk_default(color=None, msize=None, width=None):
  global _plmk_color, _plmk_width, _plmk_msize
  if color: _plmk_color = color
  if width: _plmk_width = width
  if msize: _plmk_msize = msize
  if not (color or width or msize):
    _plmk_msize = _plmk_color = _plmk_width = None

from arrayfns import *
from types import *

def spann (zmin, zmax, n = 8, fudge = 0, force = 0) :
#    return no more than N equally spaced "nice" numbers between
#    ZMIN and ZMAX.
#    HAR! I added a "force" parameter to force exactly n values
#    to be returned.
   dz = (zmax - zmin) / max (float (n), 0.)
   if dz == 0. :
      dz = abs (zmin)
   if (dz != 0.) :
      power = floor (log10 (dz) + 0.00001)
      base = dz / (10. ** power)
      if base > 5.00001 :
         base = 1.0
         power = power + 1.0
      elif base > 2.00001 :
         base = 5.0
      else :
         base = 2.0
      # Round dz up to nearest "nice" number
      dz = base * 10.0 ** power
      zmin = ceil (zmin / dz - fudge)
      zmax = floor (zmax / dz + fudge)
      if (force == 0) :
         nz = int (zmax - zmin + 1.0)
      else :
         nz = n
      if nz > 1 :
         levs = span (zmin * dz, zmax * dz, nz)
      else :
         if nz < 1 :
            if base < 1.5 :
               base = 5.0
               power = power - 1
            elif base < 2.5 :
               base = 1.0
            else :
               base = 2.0
            dz = base * 10.0 ** power
            zmin = ceil (zmin / dz + 0.001)
         levs = array ( [zmin * dz])
   else :
      levs = array ( [-1.0, 1.0])
   return (levs)

_ContourError = "ContourError"

def plfc (z, y, x, ireg, contours = 8, colors = None, region = 0,
   triangle = None, scale = "lin") :
#
# plfc (z, y, x, ireg, contours = 8, colors = None, region = 0,
#  triangle = None, scale = "lin")
#

     # 1. Get contour colors
     (vcmin, vcmax) = zmin_zmax (z, ireg)
     if type (contours) == IntType :
        n = contours
        vc = zeros (n + 2, Float)
        vc [0] = vcmin
        vc [n + 1] = vcmax
        if scale == "lin" or scale is None :
            #    This stuff is in lieu of the spann stuff in Yorick.
            vc [1:n + 1] = vcmin + arange (1, n + 1, typecode = Float) * \
               (vcmax - vcmin) / (n + 1)
        elif scale == "log" :
            vc [1:n + 1] = vcmin + exp (arange (1, n + 1, typecode = Float) * \
               log (vcmax - vcmin) / (n + 1))
        elif scale == "normal" :
            zlin = ravel (z)
            lzlin = len (zlin)
            zbar = add.reduce (zlin) / lzlin
            zs = sqrt ( (add.reduce (zlin ** 2) - lzlin * zbar ** 2) /
                (lzlin - 1))
            z1 = zbar - 2. * zs
            z2 = zbar + 2. * zs
            diff = (z2 - z1) / (n - 1)
            vc [1:n + 1] = z1 + arange (n) * diff
        else :
            raise _ContourError, "Incomprehensible scale parameter."
     elif type (contours) == ArrayType and contours.typecode () == Float :
        n = len (contours)
        vc = zeros (n + 2, Float)
        vc [0] = vcmin
        vc [n + 1] = vcmax
        vc [1:n + 1] = sort (contours)
     else :
        raise _ContourError, "Incorrect contour specification."
     if colors is None :
        colors = (arange (n + 1, typecode = Float) * (199. / n)).astype ('b')
     else :
        colors = array (colors)
        if len (colors) != n + 1 :
           raise "PLFC_Error", \
              "colors must specify one more color than contours."
        if colors.typecode != 'b' :
           colors = bytscl (colors)

     if triangle is None :
        triangle = zeros (z.shape, Int16)

     # Set mesh first
     plmesh (y, x, ireg, triangle = triangle)
     for i in range (n + 1) :
        [nc, yc, xc] = contour (array ( [vc [i], vc [i + 1]]), z)
        if (is_scalar(nc) and nc == 0 or nc is None) :
           continue
        plfp ( (ones (len (nc)) * colors [i]).astype ('b'),
           yc, xc, nc, edges = 0)
