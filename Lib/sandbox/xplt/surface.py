# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *

class Surface :
   """s= Surface ( <keyword arguments> ) will create a surface with
   coordinates and various other characteristics. The keywords are
   as follows:
        z = <value> (required). z is a two dimensional array. If x
            and y are not specified, then z will be graphed on equally
            spaced mesh points.
        x = <value>, y = <value> (if one is present, then so must
            the other be.) If c (below) is not present, this represents
            a 3d plot. Either x and y have dimensions matching z,
            or else they are one-dimensional and x's length matches
            the first dimension of z, and y's length matches the second.
        c = <value> If present, this represents a four-dimensional
            graph which is colored according to the values of c.
            c must have the same dimensions as z.
        color_card = <value> specifies which color card you wish
            to use, e. g., "rainbowhls" (the default), "random",
            etc. Although a characteristic of a Graph2d, it can
            be a Surface characteristic (of a Narcisse plot only)
            since 'link'ed surfaces can have different color cards.
        opt_3d = <value> where <value> is a string or a sequence
            of strings giving the 3d or 4d surface characteristics.
            A surface is colored by height in z if a 3d option is
            specified, and by the value of the function c if a 4d
            option is specified. With a wire grid option, the grid
            is colored; with a flat option, the quadrilaterals set
            off by grid points are colored; with a smooth option,
            the surface itself is colored by height; and with an iso
            option, the contour lines are colored. flat and iso options
            may be used together in any combination. wire grid options
            are independent of the other options. Legal arguments for
            set_3d_options are:
            'wm'--monochrome wire grid; 'w3' and 'w4'--3d and 4d
            coloring of wire grid.
            'f3' and 'f4'--flat 3d and 4d coloring options.
            'i3' and 'i4'--3d and 4d isoline (contour line) options.
            's3' and 's4'--3d and 4d smooth coloring options.
            For Gist, a subset of the above is available:
            'wm', 's3', 'f3', and 'f4'.
         ecolor = <string> In Gist, you are allowed to specify a
            color for the edges of the mesh. (default: "fg")
         shade = (0|1) Applicable only with Gist. If 1, show the
            surface shaded according to the specified or default
            light source.
         mesh_type = <string> in one of the wire modes, tells what
            form the wire grid takes: "x"--x lines only; "y": y lines only;
            "xy": both a lines and y lines.
            Not available for Gist.
         mask = <string> controls hidden line removal. Allowed values
            are "none" : transparent graph; "min": simple
            masking; "max" : better masking; "sort": slowest but
            most sophisticated.
            Gist only does "none" and "sort".
         z_c_switch = 0 or 1 : set to 1 means switch z and c in the plot.
         z_contours_scale, c_contours_scale = "lin", "log", or "normal"
         z_contours_array, c_contours_array = actual array of numbers
            to use for contours
         number_of_z_contours, number_of_c_contours = <integer>
            specifies how many contours to use; they will be computed
            automatically based on the data.
   """

   def type (self) :
      return SurfaceType

   _SurfaceSpecError = "SurfaceSpecError"

   _color_card_dict = { "absolute" : 0 , "binary" : 1 ,
      "bluegreen" : 2 , "default" : 3 , "negative" : 4 , "positive" : 5 ,
      "rainbow" : 6 , "rainbowhls" : 7 , "random" : 8 , "redblue" : 9 ,
      "redgreen" : 10 , "shifted" : 11 }

   def __init__ ( self, * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       self.cell_descr = None # (this will have values for a nonstruct. mesh)
       if ( keywords.has_key ( "z" ) ) :
          self.z = keywords ["z"]
       else :
          raise self._SurfaceSpecError , \
             "Please specify z coordinates for the surface."
       if keywords.has_key ( "x" ) :
          self.x = keywords ["x"]
          if not keywords.has_key ("y") :
             raise self._SurfaceSpecError , \
                "Both x and y coordinates must be supplied."
          self.y = keywords ["y"]
       elif keywords.has_key ( "y" ) :
          raise self._SurfaceSpecError , \
             "Both x and y coordinates must be supplied."
       else :
          self.x = None
          self.y = None
       if keywords.has_key ( "c" ) :
          self.c = keywords ["c"]
          if self.x is None : # got to dummy up x and y
             self.x = arange (self.z.shape [0], typecode = Float)
             self.y = arange (self.z.shape [1], typecode = Float)
       else :
          self.c = None
       self.generic_init ( keywords )

   def generic_init ( self , keywords ) :
       if keywords.has_key ( "color_card" ) :
          if type ( keywords [ "color_card" ] ) == type ( 1 ) :
             self.color_card = keywords [ "color_card" ]
          elif self._color_card_dict.has_key ( keywords [ "color_card" ] ) :
             self.color_card = self._color_card_dict [keywords [ "color_card" ]]
          else :
             self.color_card = 7
       else :
          self.color_card = 7
       # N. B. on the above: the color card option does not apply in Gist.
       if keywords.has_key ( "mask" ) :
          self.mask = keywords ["mask"]
       else :
          self.mask = "none"
       if keywords.has_key ( "opt_3d" ) :
          self.opt_3d = keywords ["opt_3d"]
       else :
          self.opt_3d = "wm"
       if keywords.has_key ( "mesh_type" ) :
          self.mesh_type = keywords ["mesh_type"]
       else :
          self.mesh_type = "xy"
       if keywords.has_key ( "z_c_switch" ) :
          self.z_c_switch = keywords ["z_c_switch"]
       else :
          self.z_c_switch = 0
       if keywords.has_key ( "z_contours_scale" ) :
          self.z_contours_scale = keywords ["z_contours_scale"]
       else :
          self.z_contours_scale = "lin"
       if keywords.has_key ( "c_contours_scale" ) :
          self.c_contours_scale = keywords ["c_contours_scale"]
       else :
          self.c_contours_scale = "lin"
       if keywords.has_key ( "number_of_z_contours") :
          self.number_of_z_contours = keywords ["number_of_z_contours"]
       else :
          self.number_of_z_contours = None
       if keywords.has_key ( "number_of_c_contours" ) :
          self.number_of_c_contours = keywords ["number_of_c_contours"]
       else :
          self.number_of_c_contours = None
       if keywords.has_key ( "z_contours_array" ) :
          self.z_contours_array = keywords ["z_contours_array"]
       else :
          self.z_contours_array = None
       if keywords.has_key ( "c_contours_array" ) :
          self.c_contours_array = keywords ["c_contours_array"]
       else :
          self.c_contours_array = None
       if keywords.has_key ( "ecolor" ) :
          self.ecolor = keywords ["ecolor"]
       else :
          self.ecolor = "fg"
       if keywords.has_key ( "shade" ) :
          self.shade = keywords ["shade"]
       else :
          self.shade = 0

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing surface.
       """
       del self.x, self.y, self.z, self.c, self.color_card, self.opt_3d, \
           self.mask, self.z_c_switch, self.z_contours_scale, \
           self.c_contours_scale, self.z_contours_array, \
           self.c_contours_array, self.number_of_z_contours, \
           self.number_of_c_contours, self.mesh_type
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       surface characteristics. No error checking is done.
       """
       for k in keywords.keys ():
          if k == "x" :
             self.x = keywords ["x"]
          elif k == "y" :
             self.y = keywords ["y"]
          elif k == "z" :
             self.z = keywords ["z"]
          elif k == "c" :
             self.c = keywords ["c"]
          else :
             setattr (self, k, keywords [k])
