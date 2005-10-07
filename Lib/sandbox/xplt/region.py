# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *

class Region :
   """
   r = Region ( <keyword arguments> ) is used to specify graphing modes
   for regions in a QuadMesh plot. If a QuadMesh contains one or more
   Region objects, then a Plotter, when asked to plot it, will plot
   only the regions specified.

   The keyword arguments are:

      number = <integer> the number of the Region being specified.

      boundary = 0/1 0: plot entire mesh; 1: plot only the boundary of
       the selected region. if ktype and ltype are not "none",
       then the boundary will be plotted, then the k and l lines
       with their own types.
      boundary_type, boundary_color: these matter only if boundary = 1,
       and tell how the boundary will be plotted and what its color
       will be.

      inhibit = 0/1/2 1: do not plot the (x [, j], y [, j]) lines;
        2: do not plot the (x [i, ], y[i, ]) lines.
      levels = optional two-dimensional sequences of floating point
        values. If present, a list of the values of z at which you
        want contours.
      filled = 0/1 If 1, plot a filled mesh using the values of z.
        If z is not present, the mesh zones will be filled with the
        background color, which allows plotting of a wire frame.
      edges, if nonzero when filled=1, draw a solid edge around
        each zone.
      contours = 0/1 together with filled, controls whether you get
        filled contours or just contours, or no contours at all.
      z_scale = "lin", "log", or "normal" as for QuadMesh.
      vectors = 0/1 allows the user to choose whether or not to
        plot vx and vy on this region.
      ktype, ltype: can have the same values as type, and allow the
        k and l mesh lines to be plotted differently.
      #### eventually, we can add kcolor, lcolor, kwidth, lwidth ####
      type, color, width, label, hide, marks, marker as for
        curves.
   """

   def type (self) :
       return RegionType

   _RegionSpecError = "RegionSpecError"

   def __init__ ( self , *kwds , **keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       if not keywords.has_key ( "number" ) :
          raise _RegionSpecError, "Region number not specified!"
       self.number = keywords ["number"]
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
          self.line_type = "none"
       else :
          self.line_type = "solid"
       if keywords.has_key ("type") :
          self.line_type = keywords ["type"]
       if keywords.has_key ("ktype") :
          self.ktype = keywords ["ktype"]
       else :
          self.ktype = self.line_type
       if keywords.has_key ("ltype") :
          self.ltype = keywords ["ltype"]
       else :
          self.ltype = self.line_type
       if keywords.has_key ("width") :
          self.width = keywords ["width"]
       else :
          self.width = 1
       if keywords.has_key ("color") :
          self.color = keywords ["color"]
       else :
          self.color = "fg"
       if keywords.has_key ("filled") :
          self.filled = keywords ["filled"]
       else :
          self.filled = 0
       if keywords.has_key ("edges") :
          self.edges = keywords ["edges"]
       else :
          self.edges = 0
       if keywords.has_key ("vectors") :
          self.vectors = keywords ["vectors"]
       else :
          self.vectors = 1
       if keywords.has_key ("contours") :
          self.contours = keywords ["contours"]
       elif self.filled == 0 and self.edges == 0 and self.vectors == 0 :
          self.contours = 1
       else :
          self.contours = 0
       if keywords.has_key ("z_scale") :
          self.z_scale = keywords ["z_scale"]
       else :
          self.z_scale = "lin"
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
       if keywords.has_key ("marks") :
          self.marks = keywords ["marks"]
       else :
          self.marks = 0
       if keywords.has_key ( "marker" ) :
          self.marker = keywords [ "marker" ]
       else :
          self.marker = "unspecified"

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing Region.
       """
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       Region characteristics. No error checking is done.
       """
       for k in keywords.keys ():
           if k == "type" :
              self.line_type = keywords ["type"]
           else :
              setattr (self, k, keywords [k])
