# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *
 
class Lines :
   """This class is for Gist and possibly other graphics packages;
   it creates a class of disjoint lines which can then be inserted
   into and drawn by a graph2d object.
   
   x = Lines ( <keyword arguments> ) will create a Lines object
   with coordinates and various other characteristics. The keyword
   arguments are:
       x0 = <sequence of floating point values>
       y0 = <sequence of floating point values>
       x1 = <sequence of floating point values>
       y1 = <sequence of floating point values>
          x0, y0, x1, and y1 can actually be scalars, but if arrays
          must match in size and shape. (x0[i], y0[i]) represents
          the starting point of the ith line, and (x1[i], y1[i])
          represents its endpoint.
       color = one of the legal values for Gist (currently the
          only package supporting Lines). See gist.help for details.
       hide = 0/1 (1 to hide this part of the graph)
       width = width of the lines. 1.0 (pretty narrow) is the default,
          I think successive values 2.0, 3.0, ... roughly represent
          width in pixels.
       type = "solid, "dash", "dot", "dashdot", "dashdotdot", and
          "none" (in which case the lines will be plotted as characters).
   """

   _LineSpecError = "LineSpecError"

   def type (self) :
       return LinesType

   def __init__ ( self , * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       if not keywords.has_key ("x0") or not keywords.has_key ("y0") or \
          not keywords.has_key ("x1") or not keywords.has_key ("y1" ) :
          raise self._LineSpecError , \
             "You need all four keywords x0, y0, x1, and y1."
       self.x0 = keywords ["x0"]
       self.y0 = keywords ["y0"]
       self.x1 = keywords ["x1"]
       self.y1 = keywords ["y1"]
       if len (self.x0) != len (self.y0) or len (self.x0) != len (self.x1) or \
          len (self.x0) != len (self.y1) or len (self.y0) != len (self.x1) or \
          len (self.y0) != len (self.y1) or len (self.x1) != len (self.y1) :
          raise self._LineSpecError , \
             "x0, y0, x1, and y1 must all be the same length."
       if keywords.has_key ("label" ) :
          self.label = keywords [ "label" ]
       else :
          self.label = " "
       if keywords.has_key ( "hide" ) :
          self.hide = keywords [ "hide" ]
       else :
          self.hide = 0
       if keywords.has_key ( "width" ) :
          self.width = keywords [ "width" ]
       else :
          self.width = 1.0
       if keywords.has_key ( "type" ) :
          self.line_type = keywords [ "type" ]
       else :
          self.line_type = "solid"
       if keywords.has_key ( "color" ) :
          self.color = keywords ["color"]
       else :
          self.color = "foreground" # foreground

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing curve.
       """
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       curve characteristics. No error checking is done.
       """
       for k in keywords.keys ():
          if k == "x0" :
             self.x0 = keywords ["x0"]
          elif k == "y0" :
             self.y0 = keywords ["y0"]
          elif k == "x1" :
             self.x1 = keywords ["x1"]
          elif k == "y1" :
             self.y1 = keywords ["y1"]
          elif k == "type" :
             self.line_type = keywords ["type"]
          else :
             setattr (self, k, keywords [k])
