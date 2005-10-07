# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *
from arrayfns import *
from LinearAlgebra import *
from gist import *

class QuadMesh :
   """This class is for Gist and possibly other graphics packages;
   it creates a class consisting of a QuadMesh, which can then
   be given to and plotted by a graph2d object.

   q = QuadMesh ( <keyword arguments> ) will create a QuadMesh
   object. The keyword arguments are:
   y and x, matching two-dimensional sequences of floating point
     values. These arguments are required and give the coordinates of
     the nodes of the mesh.
     --OR--
     x, y, and z can also be given as one-dimensional and
     of equal lengths. (In this case, all three arguments are
     required.) Then this will be assumed to be a random
     distribution of points (x, y) and values z on those points.
     In this case, a regular QuadMesh will be created, with nx
     equally spaced x's and ny equally spaced y's, with an nx
     by ny array of z's interpolated from the originals by
     Hardy's multiquadric fit. nx and ny can be specified by
     keyword values, lacking which they will be determined from
     ireg (if it is present), otherwise they will be defaulted
     to 50. xmax, xmin, ymax, and ymin for the new mesh may
     be specified by the corresponding keywords, or otherwise will
     by default be created by the maxima and minima of the variables
     x and y.
   ireg, optional two-dimensional sequence of integer values with
     the same dimensions as x and y, giving positive region numbers
     for the cells of the mesh, zero where the mesh does not exist.
     the first row and column of ireg are always zero, since there
     are one fewer cells in each direction than there are nodes.
   boundary = 0/1 0: plot entire mesh; 1: plot only the boundary of
     the selected region. if ktype and ltype are not "none",
     then the boundary will be plotted, then the k and l lines
     with their own types.
   boundary_type, boundary_color: these matter only if boundary = 1,
     and tell how the boundary will be plotted and what its color
     will be.
   region = n if n = 0, plot entire mesh; if any other number, plot
     the region specified (according to the settings in ireg).
   regions: "all" or a number or a list of numbers. Gives which
     regions to plot. If absent, defaults to [region], or "all"
     if region is 0. If present, and region is present also, then
     region is ignored.
   inhibit = 0/1/2 1: do not plot the (x [, j], y [, j]) lines;
     2: do not plot the (x [i, ], y[i, ]) lines.
   tri, optional two-dimensional sequence of values with
     the same dimensions as ireg, triangulation array used for
     contour plotting.
   z = optional two-dimensional sequence of floating point
     values. (Unless x and y are one-dimensional, in which case it
     must be present and match them in size. See above under x and y.)
     z has the same shape as x and y. If z is present, the
     contours of z will be plotted (default: 8 contours unless
     levels specifies otherwise), or a filled mesh will be
     plotted if filled = 1. In the latter case, z may be one
     smaller than x and y in each direction, and represents
     a zone-centered quantity.
   levels = optional one-dimensional sequence of floating point
     values. If present, a list of the values of z at which you
     want contours.
   filled = 0/1 If 1, plot a filled mesh using the values of z.
     If z is not present, the mesh zones will be filled with the
     background color, which allows plotting of a wire frame.
   contours = 0/1 only pertinent when filled = 1; if contours ia
     also 1, then will draw contours as specified by type, width,
     color, and levels.
   edges, if nonzero when filled=1, draw a solid edge around
     each zone.
   ecolor, ewidth--the color and width of the mesh lines when
     filled = 1 and adges is nonzero.
   vx, vy optional two-dimensional sequences of floating point
     values. Has the same shape as x and y. If present, represents
     a vector field to be plotted on the mesh.
   scale = floating point value. When plotting a vector field,
     a conversion factor from the units of (vx, vy) to the
     units of (x, y). 
   z_scale = specifies "log", "lin", or "normal" for how z is to be plotted.
   ktype, ltype: can have the same values as type, and allow the
     k and l mesh lines to be plotted differently.
   #### eventually, we can add kcolor, lcolor, kwidth, lwidth ####
   type, color, width, label, hide, marks, marker as for
     curves.
   """

   def type (self) :
       return QuadMeshType

   _QuadMeshSpecError = "QuadMeshSpecError"

   def __init__ ( self , * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       if not keywords.has_key ("x") or not keywords.has_key ("y") :
          raise self._QuadMeshSpecError , \
             "A QuadMesh requires both x and y keywords."
       self.x = keywords ["x"]
       self.y = keywords ["y"]
       if len (self.x.shape) == 2 :
          self.dimsofx = 2
          if self.x.shape != self.y.shape :
             raise self._QuadMeshSpecError , \
                "x and y must have the same shape."
          self.nx = self.x.shape [0]
          self.ny = self.x.shape [1]
       elif len (self.x.shape) == len (self.y.shape) == 1:
          self.dimsofx = 1
          if len (x) != len (y) :
             raise self._QuadMeshSpecError , \
                "If x and y are one dimensional, their lengths must agree."
          if keywords.has_key ("nx") :
             nx = keywords ["nx"]
          else :
             nx = None
          if keywords.has_key ("ny") :
             ny = keywords ["ny"]
          else :
             ny = None
          if keywords.has_key ("xmax") :
             xmax = keywords ["xmax"]
          else :
             xmax = max (x)
          if keywords.has_key ("ymax") :
             ymax = keywords ["ymax"]
          else :
             ymax = max (y)
          if keywords.has_key ("xmin") :
             xmin = keywords ["xmin"]
          else :
             xmin = min (x)
          if keywords.has_key ("ymin") :
             ymin = keywords ["ymin"]
          else :
             ymin = min (y)
          if keywords.has_key ("rsq") :
             rsq = keywords ["rsq"]
          else :
             rsq = 4. * (xmax - xmin) * (ymax - ymin) / len (x)
       else :
          raise self._QuadMeshSpecError , \
             "can't handle x (shape " + `self.x.shape` + ") and y " + \
             "(shape " + `self.y.shape` + "."
       if keywords.has_key ("boundary") :
          self.boundary = keywords ["boundary"]
       else :
          self.boundary = 0
       if keywords.has_key ("boundary_type") :
          self.boundary_type = keywords ["boundary_type"]
       else :
          self.boundary_type = "solid"
       if keywords.has_key ("boundary_color") :
          self.boundary_color = keywords ["boundary_color"]
       else :
          self.boundary_color = "fg"
       if keywords.has_key ("inhibit") :
          self.inhibit = keywords ["inhibit"]
       else :
          self.inhibit = 0
       if keywords.has_key ("label") :
          self.label = keywords ["label"]
       else :
          self.label = " "
       if keywords.has_key ("hide") :
          self.hide = keywords ["hide"]
       else :
          self.hide = 0
       if self.boundary == 1 :
          self.ktype = self.ltype = "none"
       else :
          self.ktype = self.ltype = "solid"
       if keywords.has_key ("type") :
          self.line_type = keywords ["type"]
       else :
          self.line_type = "solid"
       if keywords.has_key ("ktype") :
          self.ktype = keywords ["ktype"]
       if keywords.has_key ("ltype") :
          self.ltype = keywords ["ltype"]
       if keywords.has_key ("width") :
          self.width = keywords ["width"]
       else :
          self.width = 1
       if keywords.has_key ("color") :
          self.color = keywords ["color"]
       else :
          self.color = "fg"
       if keywords.has_key ("region") :
          self.region = keywords ["region"]
       else :
          self.region = 0
       if keywords.has_key ("tri") :
          self.tri = keywords ["tri"]
          if shape (tri) != shape (x) :
             raise self._QuadMeshSpecError , \
                "tri, x, and y must have the same shape."
       else :
          n1 = shape (self.x) [0]
          n2 = shape (self.x) [1]
          self.tri = zeros ( (n1, n2))
       if keywords.has_key ("ireg") and keywords ["ireg"] is not None :
          self.ireg = keywords ["ireg"]
          if self.dimsofx == 2 :
             if self.ireg.shape != self.x.shape :
                raise self._QuadMeshSpecError , \
                   "ireg, x, and y must have the same shape."
          else: # dimsofx has to be 1
             if self.nx is None :
                self.nx = self.ireg.shape [0]
             if self.ny is None :
                self.ny = self.ireg.shape [1]
             if self.ireg.shape != (nx, ny) :
                raise self._QuadMeshSpecError , \
                   "shape of ireg must be nx by ny."
       else :
          if self.dimsofx == 2:
             self.ireg = array (self.x).astype (Int)
          else :
             if self.nx is None :
                self.nx = 50
             if self.ny is None :
                self.ny = 50
             self.ireg = array ( (self.nx, self.ny), Int)
          self.ireg [0:self.nx, 0] = 0
          self.ireg [0, 0:self.ny] = 0
          self.ireg [1:self.nx, 1:self.ny] = 1
       if keywords.has_key ("filled") :
          self.filled = keywords ["filled"]
       else :
          self.filled = 0
       if keywords.has_key ("edges") :
          self.edges = keywords ["edges"]
       else :
          self.edges = 0
       if keywords.has_key ("contours") :
          self.contours = keywords ["contours"]
       elif self.filled == 0 and self.edges == 0 :
          self.contours = 1
       else :
          self.contours = 0
       if keywords.has_key ("ewidth") :
          self.ewidth = keywords ["ewidth"]
       else :
          self.ewidth = 1
       if keywords.has_key ("ecolor") :
          self.ecolor = keywords ["ecolor"]
       else :
          self.ecolor = "fg"
       if keywords.has_key ("levels") :
          self.levels = keywords ["levels"]
       else :
          self.levels = None
       if keywords.has_key ("z") :
          self.z = keywords ["z"]
       else :
          self.z = None
       if keywords.has_key ("z_scale") and keywords ["z_scale"] is not None :
          self.z_scale = keywords ["z_scale"]
       else :
          self.z_scale = "lin"
       if self.dimsofx == 1 :
          if self.z is None or len (self.z.shape) != 1 \
             or len (self.z) != len (self.x) :
             raise self._QuadMeshSpecError , \
                "If x and y are one-dimensional, " + \
                "then z must be too, and must be the same length."
          # Compute regular mesh and interpolated values
          dx = (xmax - xmin) / nx
          dy = (ymax - ymin) / ny
          xcplot = span (xmin, xmax - dx, nx)
          ycplot = span (ymin, ymax - dy, ny)
          del xmin, xmax, ymin, ymax, dx, dy
          xt = subtract.outer (self.x, self.x)
          yt = subtract.outer (self.y, self.y)
          aa = sqrt (xt * xt + yt * yt + rsq)
          del xt, yt
          alpha = array (self.z, copy = 1)
          alpha = solve_linear_equations (aa, alpha)
          self.z = mfit (alpha, self.x, xcplot, self.y,
             ycplot, rsq)
          del aa, alpha, rsq
          # Expand coordinates to 2d arrays to match zcplot
          self.x = multiply.outer (xcplot, ones (ny, Float))
          self.y = multiply.outer (ones (nx, Float), ycplot)
          del xcplot, ycplot, nx, ny
       if self.dimsofx == 2 and self.z is not None : # check for shape
          if self.filled == 0 :
             if self.z.shape != self.x.shape :
                raise self._QuadMeshSpecError , \
                   "z, x, and y must be the same shape."
          else :
             n1 = self.z.shape [0]
             n2 = self.z.shape [1]
             m1 = self.x.shape [0]
             m2 = self.x.shape [1]
             if n1 == m1 and n2 == m2 or n1 == m1 - 1 and n2 == m2 - 1 :
                pass
             else :
                raise self._QuadMeshSpecError , \
                   "z, x, and y must be the same shape, or z one smaller" + \
                   " in each dimension."
       if keywords.has_key ("vx") and not is_scalar (keywords ["vx"]) and \
          len (keywords ["vx"]) > 0 :
          if not keywords.has_key ("vy") :
             raise self._QuadMeshSpecError , \
                "vx and vy must both be present."
          self.vx = keywords ["vx"]
          self.vy = keywords ["vy"]
       else :
          self.vx = None
          self.vy = None
       if keywords.has_key ("scale") :
          self.scale = keywords ["scale"]
       else :
          self.scale = None
       if keywords.has_key ("marks") :
          self.marks = keywords ["marks"]
       else :
          self.marks = 0
       if keywords.has_key ( "marker" ) :
          self.marker = keywords [ "marker" ]
       else :
          self.marker = "unspecified"
       if keywords.has_key ( "regions" ) :
          self.regions = keywords ["regions"]
          if is_scalar (self.regions) and self.regions != "all" :
             self.regions = [self.regions]
       else :
          if self.region == 0 :
             self.regions = "all"
          else:
             self.regions = [self.region]

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing quadmesh.
       """
       self.__init__ ( keywords )
 
   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       quadmesh characteristics. Very little error checking is done.
       """
       # The following allopws vx and vy to be cleared by sending in []
       if keywords.has_key ("vx") :
          self.vx = keywords ["vx"]
          del keywords ["vx"]
       if keywords.has_key ("vy") :
          self.vy = keywords ["vy"]
          del keywords ["vy"]
       for k in keywords.keys ():
          if k == "type" :
             self.line_type = keywords ["type"]
          else :
              setattr (self, k, keywords [k])
