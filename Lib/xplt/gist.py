#  $Id$
#  ---------------------------------------------------------------------
#
#  NAME:  gist.py
#
#  Scipy CHANGES:
#  03/06/03 teo Changed all == None checks to is None checks
#  03/10/03 teo Changed eps to work with windows (noepsi=1)

#  CHANGES:
#  12/25/02 mdh Add plh to draw histograms
#  11/26/01 llc Add docstring for plmk (missing).
#  11/05/01 llc Use pydoc's help (not the one in help.py).
#  10/30/01 llc Disable is_scalar call; type(arraytype) not implemented
#               in Python.
#               Also, merge in documentation (if it exists) from gist.help 
#               after each function, so that pydoc's help can return it.
#  10/12/01 llc Re-port of gist from archived version.
#
#  ---------------------------------------------------------------------
"""
      Version: $Id$
      Copyright (c) 1996, 1997, The Regents of the University of
      California.  All rights reserved.  See Legal.htm for full
      text and disclaimer.

     Gist is a portable graphics package for scientific applications. It
     can produce interactive graphics in the X-window and Macintosh
     environments (the Python module currently supports only X11), and
     it can produce file output conforming to ANSI standard CGM or
     standard Postscript.

     Gist was developed by David H. Munro <munro@icf.llnl.gov> at
     Lawrence Livermore National Laboratory, as part of his Yorick
     scientific interpreter. Gist is distributed as part of Yorick, and
     is available at the following sites:

       wuarchive.wustl.edu: /languages/yorick/yorick-1.2.tar.gz
       sunsite.unc.edu: /pub/languages/yorick/yorick-1.2.tar.gz
       sunsite.unc.edu: /pub/Linux/apps/math/matrix/yorick-1.2.tar.gz
       netlib.att.com: /netlib/env/yorick-1.2.tar.gz
       netlib2.cs.utk.edu: /env/yorick-1.2.tar.gz

     Much of the code in the Python Gist C extension module was
     adapted from similar code in Yorick, and this help file also
     originated there. I am greatly indebted to Dave for his prior work.
     Questions about this module should be directed to me, Lee Busby,
     <busby1@llnl.gov>.

"""

__version__ = "$Id$"

from Numeric import *
from scipy_base.fastumath import *
import sys, os	# To be sure expand_path has posixpath and we have sys.path
from gistC import *
from helpmod import help as ghelp
from pydoc import help
from shapetest import *
from arrayfns import *

# Parameters used by pltitle and xytitles
pltitle_height= 18;
pltitle_font= "helvetica";

# Parameters used by plmk and plmk_default
_plmk_count = 0
_plmk_msize = _plmk_color = _plmk_width = None

#  ---------------------------------------------------------------------

def winkill(*N):
   """
   winkill( [n] )
      Delete the current graphics window, or graphics window N (0-7).
   """
   if N:
      window(N[0], display="", hcp="")
   else:
      window(display="", hcp="")

#  ---------------------------------------------------------------------

def pltitle(title):
   """
   pltitle( title )
      Plot TITLE centered above the coordinate system for any of the
      standard Gist styles.  You will need to customize this for other
      plot styles.
   """
   vp = viewport()
   xmidpt = (vp[0] + vp[1])/2.0
   plt( title, xmidpt, vp[3] + 0.02,
      font=pltitle_font, justify="CB", height=pltitle_height)

#  ---------------------------------------------------------------------

def ylimits(ymin='u',ymax='u'): 
   """
   ylimits(ymin, ymax)
      Set the y-axis plot limits in the current coordinate system to
      YMIN, YMAX, which may each be a number to fix the corresponding
      limit to a specified value, or the string "e" to make the
      corresponding limit take on the extreme value of the currently
      displayed data. Arguments may be omitted only from the right. Use
      limits( xmin, xmax ) to accomplish the same function for the x-axis
      plot limits.  Note that the corresponding Yorick function for
      ylimits is ``range'' - since this word is a Python built-in function,
      I've changed the name to avoid the collision.
      SEE ALSO: plsys, limits, logxy, plg
   """
   limits('u','u',ymin,ymax)

#  ---------------------------------------------------------------------

def moush(*arg):
   """
   moush(y, x, ireg) or moush(y, x) or moush()
      Number of args can be 0, 2, or 3.
      Returns the 1-origin zone index for the point clicked in
      for the default mesh, or for the mesh (X,Y) (region array IREG).
   """
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
      print "Mouse takes 0, 2, or 3 args: ( [ y, x [ , ireg ] ] )"
      return None

#  ---------------------------------------------------------------------
#  PURPOSE:  Create an encapsulated PostScript file.  
#            Requires Ghostscript and its associated ps2epsi utility.
#  ---------------------------------------------------------------------

def eps(name, epsi=0, pdf=0):
   """
      Write the picture in the current graphics window to the Encapsulated
      PostScript file NAME+".eps" (i.e.- the suffix .eps is added to NAME).
      The last extension of name is stripped to avoid .eps.eps files
      
      If epsi is 1, this function requires the ps2epsi utility which
      comes with the project GNU Ghostscript program and will place a
      bitmap image in the postscript file. 
      Any hardcopy file associated with the current window is first closed,
      but the default hardcopy file is unaffected.  As a side effect,
      legends are turned off and color table dumping is turned on for
      the current window.
      SEE ALSO: window, fma, hcp, hcp_finish, plg
   """
   import os
   name,ignore = os.path.splitext(name)
   totalname = name
   apath, basename = os.path.split(name)
   name = name + ".ps"   
   window (hcp = name, dump = 1, legends = 0)
   hcp ()
   window (hcp="")
   res = 1
   if epsi:
      res = os.system ("ps2epsi " + name)
   if not res:
      os.remove(name)
      os.rename ("%s.epsi" % basename, "%s.eps" % totalname)
   else:
      os.rename("%s.ps" % totalname, "%s.eps" % totalname)
   if pdf:
      os.system("epstopdf %s.eps" % totalname)

#  ---------------------------------------------------------------------

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

#  ---------------------------------------------------------------------

def _spanz(lb,ub,n):

   if n < 3: raise ValueError, '3rd arg must be at least 3'
   c = 0.5*(ub - lb)/(n - 1.0)
   b = lb + c
   a = (ub - c - b)/(n - 2.0)
   return map(lambda x,A=a,B=b: A*x + B, arange(n-1))

#  .. predefined markers: square, +, delta, circle, diamond, x, grad

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

#  ---------------------------------------------------------------------

def plmk(y,x=None,marker=None,width=None,color=None,msize=None):
   """
   Make a scatter plot of the points Y versus X.  If X is nil,
   it defaults to indgen(numberof(Y)).  By default, the marker
   cycles through 7 predefined marker shapes.  You may specify a shape
   using the marker= keyword, line width using the width= keyword (you
   get solid fills for width>=10), color using the color= keyword.
   You can also use the msize= keyword to scale the marker (default
   msize=1.0).  You can change the default width, color, or msize
   using the plmk_default function.

   The predefined marker= values are:

   marker=
      1        square
      2        cross
      3        triangle
      4        circle
      5        diamond
      6        cross (rotated 45 degrees)
      7        triangle (upside down)

   You may also put marker=[xm,ym] where xm and ym are vectors
   of NDC coordinates to design your own custom marker shapes.

   SEE ALSO: plmk_default, plg (type=0 keyword)
   """

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

#  ---------------------------------------------------------------------


def plmk_default(color=None, msize=None, width=None):
   """
   Set default color, msize, and width values for plmk.  Use
   width=10 to get solid fills.  With no parameters, plmk_default
   restores the initial default values.
   """
   global _plmk_color, _plmk_width, _plmk_msize
   if color: _plmk_color = color
   if width: _plmk_width = width
   if msize: _plmk_msize = msize
   if not (color or width or msize):
      _plmk_msize = _plmk_color = _plmk_width = None

from arrayfns import *
from types import *

#  ---------------------------------------------------------------------

def spann (zmin, zmax, n = 8, fudge = 0, force = 0) :
   """
   spann (zmin, zmax, n = 8, fudge = 0, force = 0)
      return no more than N equally spaced "nice" numbers between
      ZMIN and ZMAX.
      Note that in general spann may not supply the number of
      values that you asked for. To force it to do so, set
      keyword FORCE to nonzero.
      SEE ALSO: span, spanl, plc, plfc
   """
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

#  ---------------------------------------------------------------------

def plfc (z, y, x, ireg, contours = 8, colors = None, region = 0,
   triangle = None, scale = "lin") :
   """
   plfc (z, y, x, ireg, contours = 8, colors = None, region = 0,
      triangle = None, scale = "lin")
      fills contours of Z on the mesh Y versus X.  Y, X, and IREG are
      as for plm.  The Z array must have the same shape as Y and X.
      The function being contoured takes the value Z at each point
      (X, Y) -- that is, the Z array is presumed to be point-centered.

      The CONTOURS keyword can be an integer specifying the number of
      contours desired, or a list of the values of Z at which you want
      contour curves.  These curves divide the mesh into len(CONTOURS+1)
      regions, each of which is filled with a solid color.  If CONTOURS is
      None or not given, 8 "nice" equally spaced level values spanning the
      range of Z are selected.

      If you specify CONTOURS, you may also specify COLORS, an array of
      color numbers (Python typecode 'b', integers between 0 and the
      length of the current palette - 1, normally 199) of length
      len(CONTOURS)+1. If you do not specify them, equally
      spaced colors are chosen.

      If CONTOURS is an integer, SCALE expresses how contour levels
      are determined.  SCALE may be "lin", "log", or "normal"
      specifying linearly, logarithmically, or normally spaced
      contours. Note that unlike Yorick's plfc, this routine does
      not use spann to compute its contours. Neither, apparently,
      does plc, which uses a third algorithm which matches neither
      the one we use nor the one spann uses. So if you plot filled
      contours and then plot contour lines, the contours will in
      general not coincide exactly.

      Note that you may use spann to calculate your contour levels
      if you wish.

      The following keywords are legal (each has a separate help entry):
    KEYWORDS: triangle, region
    SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh
              color_bar, spann, contour, limits, logxy, range, fma, hcp
   """
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
#     is_scalar and type  do not work with multiarray objects or rank-0 arrays
#     .. Try skipping
#     if (is_scalar(nc) and nc == 0 or nc is None) :
#        continue
      plfp ( (ones (len (nc)) * colors [i]).astype ('b'),
         yc, xc, nc, edges = 0)

def plh (y, x = None, width = 1, hide = 0, legend = None, color = None):
   """
   plh ( y, x = None)
      draws a histogram, where the height of the bars is given
      by Y. If X is None, the bars in the histogram have a
      width equal to unity. If X is a single real number, the
      widths of the bars are all equal to this number. If X is
      a one-dimensional array with an equal number of elements
      as the Y array, then X is interpreted as the widths of
      the bars. In all of these cases, the histogram starts at
      the origin. However, if X is a one-dimensional array
      with one element more than Y, then X is interpreted as
      the locations of the start and end points of the bars in
      the histogram.

      The following keywords are legal (each has a separate help entry):
    KEYWORDS: legend, width, hide, color
    SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmesh
              color_bar, spann, contour, limits, logxy, range, fma, hcp
   """

   n = len(y)
   barx = zeros(n*4,'d')
   if x==None:
      for i in range(n-1):
         barx[4*i+2:4*i+6] = i+1
      barx[-2:] = n
   else:
      if type(x) == IntType or type(x) == FloatType:
         for i in range(n-1):
            barx[4*i+2:4*i+6] = i * x
         barx[-2:] = n * x
      elif type(x) == ListType or type(x) == ArrayType:
         if len(x) == n:
            for i in range(n-1):
               barx[4*i+2:4*i+6] = barx[4*i] + x[i]
            barx[-2:] = barx[-3] + x[-1]
         elif len(x) == n + 1:
            barx[:2] = x[0]
            for i in range(1,n+1):
               barx[4*i-2:4*i+2] = x[i]
            barx[-2:] = x[-1]
         else:
            raise "plh error: inconsistent length of X"
   bary = zeros(4*n,'d')
   for i in range(n):
      bary[4*i+1:4*i+3] = y[i]
   pldj (barx[1:], bary[1:], barx[:-1], bary[:-1], width=width, hide=hide, \
         legend = legend, color = color)

