# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *
from yorick import *
 
class Curve :

   """
   x = Curve (...keyword arguments...) will create a curve with coordinates
   and characteristics such as color, type, label, etc.
   The keyword arguments for Curve are:
 
       y = <sequence of floating point values> (required)
       x = <sequence of floating point values> (optional)
       color = <value> where <value> is an integer from 0 to 63
               representing an entry in a color chart, or a
               common color name like "red", "blue", "background", etc.
       axis = "left" or "right" tells whether the left or right
              y axis will be assigned to this curve. (This option is
              not available in all graphics packages and will be ignored
              if that is the case.)
       label = <string> represents the label of this curve.
       type = <value> tells how the curve will be plotted: "line",
              "solid" (same as "line"), "step", "dash", "dashdot".
              "dashdotdot", "none", "+", "*", "o", "x", and "." are
              allowed. If the option is not available in a particular
              graphics package, a good guess will be substituted.
              If type = "none" and marks = 1, the plot will be
              a polymarker plot, if supported by the graphics.
              Note that because of disparities among graphics packages
              supported, you can specify plotting a curve
              pointwise with symbols like "+", "*", etc., either
              by use of the type variable or using marks and markers
              in conjunction with type = "none".
       marks = 0 or 1; select unadorned lines (0) or lines with
              occasional markers (1). Some graphics packages do not
              have this option. The markers default to letters of
              the alphabet, but can be changed by the marker keyword.
       marker = character or integer value for character used to mark
              this curve if marks = 1. Special values '\1', '\2',
              '\3', '\4', and '\5' stand for point, plus, asterisk,
              circle, and cross, which sometimes look prettier than
              characters on some devices.
       width = real number; specifies the width of a curve if this
              is supported by the graphics. 1.0 gives a pretty finely
              drawn curve and is the default.
       hide = 0 or 1; if set to 1, this curve will be hidden on
              the plot.
 
   Let x be a Curve object. Then some of the methods are:
 
   x.set (...keyword arguments...) will set the specified argument(s)
      to new value(s)
   x.new (...keyword arguments...) clears the curve and defines a
      new one.
   """

   def type (self) :
       return CurveType

   _CurveSpecError = "CurveSpecError"

   def __init__ ( self, * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       if ( keywords.has_key ( "y" ) ) :
          self.y = keywords ["y"]
       else :
          raise self._CurveSpecError , "No y array has been specified."
       if is_scalar (self.y) :
          raise self._CurveSpecError , "y must be an array."
       if ( keywords.has_key ( "x" ) and keywords ["x"] is not None) :
          self.x = keywords ["x"] 
       else :
          if len ( shape (self.y) ) == 1 :
             self.x = arange (len (self.y), typecode = Float)
          else :
             self.x = arange (self.y.shape [1], typecode = Float)
       if abs (len (self.x) - len (self.y)) == 1 :
          # average the longer one
          if len (self.x) > len (self.y) :
             self.x = zcen_ (self.x, 0)
          else :
             self.y = zcen_ (self.y, 0)
       elif len (self.x) != len (self.y) :
          raise self._CurveSpecError , \
             "x and y must satisfy |len (x) - len (y)| <= 1."
       if keywords.has_key ( "type" ) :
          # ZCM 5/23/96 adding gist, pass type on unchanged
          self.line_type = keywords [ "type" ]
       else :
          self.line_type = "unspecified"
       # New type added for Gist
       if keywords.has_key ( "marker" ) :
          self.marker = keywords [ "marker" ]
       else :
          self.marker = None # default on engines with marked curves
       if keywords.has_key ( "marks" ) :
          self.marks = keywords [ "marks" ]
       else :
          self.marks = 0
       if keywords.has_key ( "color" ) :
          # ZCM 5/23/96 color conversions done at curve2d level
          self.color = keywords ["color"]
       else :
          self.color = "foreground" # foreground
       if keywords.has_key ( "axis" ) :
          self.axis = keywords [ "axis" ]
       else :
          self.axis = "left"
       if keywords.has_key ( "label" ) :
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
       if keywords.has_key ( "marker" ) :
          self.marker = keywords [ "marker" ]
       else :
          self.marker = None # default on engines with marked curves
       if keywords.has_key ( "marks" ) :
          self.marks = keywords [ "marks" ]
       else :
          self.marks = 0

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing curve.
       """
       del self.x, self.y, self.axis, self.label, self.line_type, self.color
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       curve characteristics. No error checking is done.
       """
       for k in keywords.keys ():
          if k == "x" :
             self.x = keywords ["x"]
          elif k == "y" :
             self.y = keywords ["y"]
          elif k == "type" :
             self.line_type = keywords ["type"]
          else :
             setattr (self, k, keywords [k])
       if abs (len (self.x) - len (self.y)) == 1 :
          # average the longer one
          if len (self.x) > len (self.y) :
             self.x = zcen_ (self.x, 0)
          else :
             self.y = zcen_ (self.y, 0)
       elif len (self.x) != len (self.y) :
          raise self._CurveSpecError , \
             "x and y must satisfy |len (x) - len (y)| <= 1."
