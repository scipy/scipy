# $Id$
# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

#
#  PLWF.PY
#  Simple "painter's algorithm"-class routine for making 3-D wire frames
#  and related models.
#
#  $Id$
#

## execfile ("pl3d.py")
from types import *
from arrayfns import *
from pl3d import *

def plwf (z, y = None, x = None, fill = None, shade = 0, edges = 1,
   ecolor =  None, ewidth = None, cull = None, scale = None, cmax = None,
   clear = 1) :

   """
   plwf (z)
   or plwf (z, y, x)

     plots a 3-D wire frame of the given Z array, which must have the
     same dimensions as the mesh (X, Y).  If X and Y are not given, they
     default to the first and second indices of Z, respectively.
     The drawing order of the zones is determined by a simple "painter's
     algorithm", which works fairly well if the mesh is reasonably near
     rectilinear, but can fail even then if the viewpoint is chosen to
     produce extreme fisheye perspective effects.  Look at the resulting
     plot carefully to be sure the algorithm has correctly rendered the
     model in each case.

   KEYWORDS: fill   -- optional colors to use (default is to make zones
                       have background color), same dimension options as
                       for z argument to plf function
             shade  -- set non-zero to compute shading from current
                       3D lighting sources
             edges  -- default is 1 (draw edges), but if you provide fill
                       colors, you may set to 0 to supress the edges
             ecolor, ewidth  -- color and width of edges
             cull   -- default is 1 (cull back surfaces), but if you want
                       to see the "underside" of the model, set to 0
             scale  -- by default, Z is scaled to "reasonable" maximum
                       and minimum values related to the scale of (X,Y).
                       This keyword alters the default scaling factor, in
                       the sense that scale=2.0 will produce twice the
                       Z-relief of the default scale=1.0.
             cmax   -- the ambient= keyword in light3 can be used to
                       control how dark the darkest surface is; use this
                       to control how light the lightest surface is
                       the lightwf routine can change this parameter
                       interactively

   SEE ALSO: lightwf, plm, plf, orient3, light3, fma3, window3
   """

   _draw3 = get_draw3_ ( )
   _square = get_square_ ( )
   [_xfactor, _yfactor] = get_factors_ ( )

   if (type (z) == ListType) :
      xyz = z [0]
      fill = z [1]
      shade = z [2]
      edges = z [3]
      ecolor = z [4]
      ewidth = z [5]
      cull = z [6]
      cmax = z [7]
      
      xyz1 = get3_xy(xyz, 1)
      x = xyz [0] # the original x
      y = xyz [1] # the original y
      

      # rotate (x,y,0) into on-screen orientation to determine order
      # just use four corners for this
      nx = shape (x)
      ny = nx [1]
      nx = nx [0]
      xx = array([[x [0, 0], x[nx - 1, 0]],
                  [x [0, ny - 1] , x[nx - 1, ny - 1]]])
      yy = array([[y [0, 0], y[nx - 1, 0]],
                  [y [0, ny - 1] , y[nx - 1, ny - 1]]])
      xyzc = array ( [ xx , yy, array ( [ [0., 0.], [0., 0.]])])
      xyzc = get3_xy(xyzc, 1)

      # compute mean i-edge and j-edge vector z-components
      iedge = avg_ (xyzc [2, :, -1] - xyzc [2, :, 0])
      jedge = avg_ (xyzc [2, -1] - xyzc [2, 0])

      # compute shading if necessary
      if (shade) :
         xyz = xyz1
         fill = get3_light (xyz)
      # The order either requires a transpose or not, reversal of the
      # order of the first dimension or not, and reversal of the order
      # of the second dimension or not.

      # The direction with the minimum magnitude average z-component must
      # vary fastest in the painting order.  If this is the j-direction,
      # a transpose will be required to make this the i-direction.
      if abs (iedge) < abs (jedge) :
         tmp = iedge
         iedge = jedge
         jedge = tmp
         x = transpose (array (xyz1 [0]))
         y = transpose (array (xyz1 [1]))
         if fill != None :
            fill = transpose (fill)
      else :
         x = xyz1 [0]
         y = xyz1 [1]

      # Zones must be drawn from back to front, which means that the
      # average z-component of the edge vectors must be positive.  This
      # can be arranged by reversing the order of the elements if
      # necessary.
      if iedge < 0.0 :
         x = reverse (x, 0)
         y = reverse (y, 0)
         if fill != None :
            fill = reverse (fill, 0)
      if jedge < 0.0 :
         x = reverse (x, 1)
         y = reverse (y, 1)
         if fill != None :
            fill = reverse (fill, 1)
      xmax = maxelt_ (x)
      xmin = minelt_ (x)
      ymax = maxelt_ (y)
      ymin = minelt_ (y)
      if _xfactor != 1. :
         xmax = xmax + (_xfactor - 1) * (xmax - xmin) / 2.0
         xmin = xmin - (_xfactor - 1) * (xmax - xmin) / 2.0
      if _yfactor != 1. :
         ymax = ymax + (_yfactor - 1) * (ymax - ymin) / 2.0
         ymin = ymin - (_yfactor - 1) * (ymax - ymin) / 2.0
      if _square :
         xdif = xmax - xmin
         ydif = ymax - ymin
         if xdif > ydif :
            dif = (xdif - ydif) / 2.
            ymin = ymin - dif
            ymax = ymax + dif
         elif ydif > xdif :
            dif = (ydif - xdif) / 2.
            xmin = xmin - dif
            xmax = xmax + dif
      if fill != None :
         if len (fill.shape) == 1:
            fill = bytscl (fill)
         else:
            k = fill.shape [0]
            l = fill.shape [1]
            fill = reshape ( bytscl (ravel (fill)), (k, l))
      if cull == 0 : #transparent mesh
         if ecolor != None :
            plm (y, x, color = ecolor)
         else :
            plm (y, x)
      elif ecolor != None and ewidth != None and cmax != None :
         plf (fill, y, x, edges = edges, ecolor = ecolor,
              ewidth = ewidth, cmin = 0.0, cmax = cmax, legend = "")
      elif ecolor != None and ewidth != None :
         plf (fill, y, x, edges = edges, ewidth = ewidth,
              cmin = 0.0, ecolor = ecolor, legend = "")
      elif ecolor != None and cmax != None :
         plf (fill, y, x, edges = edges, ecolor = ecolor,
              cmin = 0.0, cmax = cmax, legend = "")
      elif ewidth != None and cmax != None :
         plf (fill, y, x, edges = edges,  ewidth = ewidth,
              cmin = 0.0, cmax = cmax, legend = "")
      elif ecolor != None :
         plf (fill, y, x, edges = edges, ecolor = ecolor,
              cmin = 0.0, legend = "")
      elif ewidth != None :
         plf (fill, y, x, edges = edges, ewidth = ewidth,
              cmin = 0.0, legend = "")
      elif cmax != None :
         plf (fill, y, x, edges = edges,
              cmin = 0.0, cmax = cmax, legend = "")
      else :
         plf (fill, y, x, edges = edges, cmin = 0.0, legend = "")
      return [xmin, xmax, ymin, ymax]

   xyz = xyz_wf (z, y, x, scale = scale)

   if clear :
      clear3 ( )
   set3_object (plwf, [xyz, fill, shade, edges, ecolor, ewidth, cull, cmax])
   if ( _draw3 ) :
      call_idler ( ) # This will traverse and execute the drawing list
                     # if the default idler has been set.

_LightwfError = "LightwfError"

def lightwf (cmax) :

   """
   lightwf (cmax)
     Sets the cmax= parameter interactively, assuming the current
     3D display list contains the result of a previous plwf call.
     This changes the color of the brightest surface in the picture.
     The darkest surface color can be controlled using the ambient=
     keyword to light3.

   SEE ALSO: plwf, light3
   """

   _draw3_list = get_draw3_list_ ()
   _draw3_n = get_draw3_n_ ()
   list = _draw3_list [_draw3_n:]
   if list [0] != plwf :
      raise _LightwfError, "current 3D display list is not a plwf"
   list [1] [7] = cmax
   undo3_set_ (lightwf, list)


_Xyz_wfError = "Xyz_wfError"

def xyz_wf (z, y, x, scale = 1.0) :

   """
   xyz_wf (z, [y, x] [,scale = 1.0])
      returns a 3-by-ni-by-nj array whose 0th entry is x, 1th entry
      is y, and 2th entry is z. z is ni-by-nj. x and y, if present,
      must be the same shape. If not present, integer ranges will
      be used to create an equally spaced coordinate grid in x and y.
      The function which scales the "topography" of z(x,y) is
      potentially useful apart from plwf.
      For example, the xyz array used by plwf can be converted from
      a quadrilateral mesh plotted using plf to a polygon list plotted
      using plfp like this:
        xyz= xyz_wf(z,y,x,scale=scale);
        ni= shape(z)[1];
        nj= shape(z)[2];
        list = ravel (add.outer (
           ravel(add.outer (adders,zeros(nj-1, Int))) +
           arange((ni-1)*(nj-1), typecode = Int),
           array ( [[0, 1], [nj + 1, nj]])))
        xyz=array([take(ravel(xyz[0]),list),
           take(ravel(xyz[1]),list),
           take(ravel(xyz[2]),list)])
        nxyz= ones((ni-1)*(nj-1)) * 4;
      The resulting array xyz is 3-by-(4*(nj-1)*(ni-1)).
      xyz[0:3,4*i:4*(i+1)] are the clockwise coordinates of the
      vertices of cell number i.
   """

   if len (shape (z)) < 2 :
      raise _Xyz_wfError, "impossible dimensions for z array"
   nx = shape (z) [0]
   ny = shape (z) [1]
   if y == None or x == None :
      if x != None or y != None :
         raise _Xyz_wfError, "either give y,x both or neither"
      x = span (0, ny - 1, ny, nx)
      y = transpose (span (0, nx - 1, nx, ny))
   elif shape (x) != shape (z) or shape (y) != shape (z) :
      raise _Xyz_wfError, "x, y, and z must all have same dimensions"
   xyscl = max (maxelt_ (x) - minelt_ (x),
                maxelt_ (y) - minelt_ (y))
   if scale != None:
      xyscl = xyscl * scale
   dz = maxelt_ (z) - minelt_ (z)
   zscl= dz + (dz == 0.0)
   if zscl :
      z = z * 0.5 * xyscl /zscl
   xbar = avg_ (x)
   ybar = avg_ (y)
   zbar = avg_ (z)
   xyz = array ( [x - xbar, y - ybar, z - zbar], Float)
   return (xyz)
