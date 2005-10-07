# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from surface import *
from graftypes import *
from plane import *
from slice3 import *
 
class Mesh3d ( Surface ) :
   """m = Mesh3d ( <keyword arguments> ) will create a mesh with
   coordinates and various other characteristics. The keywords are
   as follows:
      For all meshes:
         color_card = <value> specifies which color card you wish
            to use, e. g., "rainbowhls" (the default), "random",
            etc. Although a characteristic of a Graph2d, it can
            be a Surface characteristic since 'link'ed surfaces
            can have different color cards.
            This is not implemented for Gist. Instead, Gist
            allows the option of a "split palette", which causes
            isosurfaces to be shaded and other surfaces to be
            colored by a subset of the current (or specified)
            palette.
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
            For Gist plots, currently only 'wm', 'f3', and 'f4'
            are currently available. If others are specified, an
            intelligent guess will be made.
         mesh_type = <string> in one of the wire modes, tells what
            form the wire grid takes: "x"--x lines only; "y": y lines only;
            "xy": both a lines and y lines.
            Not implemented for Gist plots.
         mask = <string> controls hidden line removal. Allowed values
            are "none" : transparent graph; "min": simple
            masking; "max" : better masking; "sort": slowest but
            most sophisticated.
            Only the "none" and "sort" are available for Gist.
         z_c_switch = 0 or 1 : set to 1 means switch z and c in the plot.
         z_contours_scale, c_contours_scale = "lin" or "log"
         z_contours_array, c_contours_array = actual array of numbers
            to use for contours
         number_of_z_contours, number_of_c_contours = <integer>
            specifies how many contours to use; they will be computed
            automatically based on the data.
         The last six options are not available in Gist.
      For nonstructured meshes only:
         x = <values> , y = <values> , z = <values> three vectors
            of equal lengths giving the coordinates of the nodes of
            a nonstructured mesh.
         c = <values> , a vector of the same size as x, y, and z
            giving a data value at each of the points specified.
            c could also be an array of such vectors, when
            isosurfaces of more than one function are to be plotted.
         avs = 0 or 1: if 1, the input data represents a nonstructured
            mesh in a sort of AVS format, which will be explained in
            more detail below. The data will be translated into the
            Narcisse format prior to being sent to Narcisse.
         cell_descr = <integer array> if present, this keyword signifies
            a nonstructured mesh submitted in the Narcisse format,
            as explained in the Narcisse manual. For nonstructured
            mesh data, one of these two keywords must be present.
         In the case avs = 1, at least one of the following four keywords
            (hex, tet, prism, or pyr) must be present also:
         hex = [ <list of hexahedral cell data> ] The entries in the
            list (in order) must be:
            an integer number n_zones, which is the number of hex
            cells in the mesh; and a matrix nz whose dimensions are
            n_zones by 8; nz [i][0], nz [i][1], ... nz [i][7] give
            the indices of the 8 vertices of the ith zone in canonical
            order (one side in the outward normal direction, then the
            corresponding vertices of the opposite side in the inward
            normal direction).
         tet = [ <list of tetrahedral cell data> ] The list is the
            same format as for hex data. The matrix nz will be
            n_zones by 4, and each row gives the indices of the
            apex and then the base in inward normal order.
         prism = [ <list of prismatic cell data> ] The list is the
            same format as for hex data. The matrix nz will be
            n_zones by 6, and each row gives the indices of one
            triangular side in the outward normal direction, then
            the corresponding vertices of the opposite side in the
            inward normal direction.
         pyr = [ <list of pyramidal cell data> ] The list is the
            same format as for hex data. The matrix nz will be
            n_zones by 5, and each row gives the indices of the
            apex and then the base in inward normal order.
         In the case of avs = 0, the following keywords are all
            required:
         cell_descr = <values> , a vector of integers giving the
            description of the cells of the mesh, as specified in
            the Narcisse manual.
         no_cells = <integer value> , the number of cells in the mesh.
            (Here cells refers to the number of faces of cells).
      For structured meshes only:
         x = <values> , y = <values> , z = <values> three vectors
            giving the coordinates of the nodes of a structured
            rectangular mesh.
         c = <values> , a three-dimensional array whose [i, j, k]
            element gives the associated data value at (x [i], y [j],
            z [k]). Could also be one less in each direction,
            giving a cell-centered quantity.
            c could also be an array of such arrays, when
            isosurfaces of more than one function are to be plotted.
   """

   _MeshInitError = "MeshInitError"

   def type (self) :
       return Mesh3dType

   def __init__ ( self , *kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       self.cell_dict = {}
       if keywords.has_key ("hex") :
          self.cell_dict ["hex"] = keywords ["hex"]
       if keywords.has_key ("tet") :
          self.cell_dict ["tet"] = keywords ["tet"]
       if keywords.has_key ("prism") :
          self.cell_dict ["prism"] = keywords ["prism"]
       if keywords.has_key ("pyr") :
          self.cell_dict ["pyr"] = keywords ["pyr"]
       if (keywords.has_key ("avs") and keywords ["avs"] == 0) or \
          (not keywords.has_key ("avs")
           and not keywords.has_key ("cell_descr")) :
          # structured case works basically just like a Surface
          Surface.__init__ ( self, keywords )
          self.structured = 1
       else : # nonstructured case, at least initialize the generics
          self.structured = 0
          self.generic_init ( keywords )
          if keywords.has_key ("avs") :
             self.avs = keywords ["avs"]
          else :
             self.avs = None
          if not keywords.has_key ("x") or \
             not keywords.has_key ("y") or \
             not keywords.has_key ("z") or \
             not keywords.has_key ("c") :
             raise self._MeshInitError , \
                "x, y, z, and c keywords required for nonstructured mesh."
          self.x = keywords ["x"]
          self.y = keywords ["y"]
          self.z = keywords ["z"]
          self.c = keywords ["c"]
          if keywords.has_key ("cell_descr" ) : # already there
             if not keywords.has_key ("no_cells") :
                raise self._MeshInitError , \
                   "Keyword no_cells is required if avs is not specified."
             self.cell_descr = keywords ["cell_descr"]
             self.number_of_cells = keywords ["no_cells"]
          elif not keywords.has_key ("avs") :
             raise self._MeshInitError , \
                "keywords must include avs = 1 or else cell_descr = values."

   # Conversion tables for AVS to Narcisse formats

   _cell_types = ["tet", "pyr", "prism", "hex"]
   _no_of_vertices = { "hex" : 8, "tet" : 4, "pyr" : 5, "prism" : 6 }
   _no_of_rect_faces = { "hex" : 6, "tet" : 0, "pyr" : 1, "prism" : 3 }
   _no_of_tri_faces = { "hex" : 0, "tet" : 4, "pyr" : 4, "prism" : 2 }
   # The next few items give the standard numbering of the faces of
   # the various cells.
   _hex_triangular_faces = [ ]
   _tet_triangular_faces = [ [0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]]
   _pyr_triangular_faces = [ [0, 1, 2], [0, 4, 1], [0, 3, 4], [0, 2, 3]]
   _prism_triangular_faces = [ [0, 1, 2], [3, 5, 4]]
   _hex_rectangular_faces = [ [0, 1, 2, 3], [1, 0, 4, 5], [2, 1, 5, 6],
                              [3, 2, 6, 7], [0, 3, 7, 4], [4, 7, 6, 5]]
   _tet_rectangular_faces = [ ]
   _pyr_rectangular_faces = [ [1, 4, 3, 2]]
   _prism_rectangular_faces = [ [0, 3, 4, 1], [1, 4, 5, 2], [2, 5, 3, 0]]
   _tri_face_numbers = { "hex" : _hex_triangular_faces , 
                         "tet" : _tet_triangular_faces ,
                         "pyr" : _pyr_triangular_faces ,
                         "prism" : _prism_triangular_faces }
   _rect_face_numbers = { "hex" : _hex_rectangular_faces ,
                          "tet" : _tet_rectangular_faces ,
                          "pyr" : _pyr_rectangular_faces ,
                          "prism" : _prism_rectangular_faces }

   def get_verts_list (self) :
      """get_verts_list () returns a list of vertex arrays in the
      form wanted by mesh3. Note: if they are in avs order,
      they need to be put into Gist order (pyramids and tets
      are OK, only prisms and hexahedra need vertices permuted).
      """
      verts = []
      for k in self._cell_types :
         if self.cell_dict.has_key (k) :
            n_z = self.cell_dict[k][1]
            if n_z.shape [1] == 8 :
               verts.append (
                  take (n_z, array ( [0, 1, 3, 2, 4, 5, 7, 6])))
            elif n_z.shape [1] == 6 :
               verts.append (
                  take (n_z, array ( [3, 0, 4, 1, 5, 2])))
            else :
               verts.append (n_z)
      if len (verts) == 1 :
         verts = verts [0]
      return verts

   def create_Narcisse_format ( self ) :
      """create_Narcisse_format ( <keyword dict> ) is an internal
      routine which takes input in (essentially) avs format and
      converts it into the format desired by Narcisse.
      """
      if hasattr (self, "cell_descr") : # don't calculate it again.
         return
      if not self.cell_dict.has_key ("hex") \
         and not self.cell_dict.has_key ("tet") \
         and not self.cell_dict.has_key ("prism") \
         and not self.cell_dict.has_key ("pyr") :
         raise _MeshInitError, \
            "In avs mode, at least one of the keywords " + \
            "(hex, tet, pyr, prism) is required."
      n_zones = { }
      nz = { }
      for k in self._cell_types :
         if self.cell_dict.has_key (k) :
            n_zones [k] = self.cell_dict [k][0]
            nz [k] = self.cell_dict[k][1]
         else :
            n_zones [k] = 0
            nz [k] = [ ]
      self.cell_descr = zeros (30 * n_zones ["hex"] +
                               21 * n_zones ["pyr"] +
                               16 * n_zones ["tet"] +
                               23 * n_zones ["prism"], Int)
      self.number_of_cells = n_zones ["hex"] * 6 + n_zones ["pyr"] * 5 + \
                             n_zones ["tet"] * 4 + n_zones ["prism"] * 5
      n = 0
      start = 0
      # Fill the start of the cell_descriptor array
      # The first part is simply a sequence of indices into the second
      # part of the self.cell_descr array; Each entry tells where the indices
      # for this particular face begin.
      cell_j = 0
      for k in self._cell_types :
         nrf = self._no_of_rect_faces [k]
         ntf = self._no_of_tri_faces [k]
         for i in range (start, start + n_zones [k]) :
            for j in range (ntf) : # count triangular faces first
               n = n + 3
               self.cell_descr [cell_j] = n
               cell_j = cell_j + 1
            for j in range (nrf) : # count rectangles next
               n = n + 4
               self.cell_descr [cell_j] = n
               cell_j = cell_j + 1
         start = start + n_zones [k]
      # Next we go through each type of cell and enter the indices
      # of its faces, face-by-face, triangular sides first by convention.
      for k in self._cell_types :
         nrf = self._no_of_rect_faces [k]
         ntf = self._no_of_tri_faces [k]
         n_z = nz [k]
         for i in range (start, start + n_zones [k]) :
            n_z_ = n_z [i - start]
            for j in range (ntf) : # do triangular faces first
               ind = self._tri_face_numbers [k][j]
               for l in range (3) :
                  self.cell_descr [cell_j + l] = n_z_ [ind [l]]
               cell_j = cell_j + 3
            for j in range (nrf) : # do rectangular faces next
               ind = self._rect_face_numbers [k][j]
               for l in range (4) :
                  self.cell_descr [cell_j + l] = n_z_ [ind [l]]
               cell_j = cell_j + 4
         start = start + n_zones [k]

   def new ( self , ** keywords ) :
      """new (...keyword arguments...) allows you to reuse a
      previously defined mesh.
      """
      del self.x, self.y, self.z, self.c, self.color_card, self.opt_3d, \
          self.mask, self.z_c_switch, self.z_contours_scale, \
          self.c_contours_scale, self.z_contours_array, \
          self.c_contours_array, self.number_of_z_contours, \
          self.number_of_c_contours, self.mesh_type
      if not self.structured :
         del self.cell_descr, self.number_of_cells
      self.__init__ ( keywords )

from plane import *

class Slice :
   """
   Class Slice is used to contain a geometric description of
   an isosurface or plane Slice of a Mesh3d object. It is 
   created by the builtin function sslice, defined later in
   this module.
   Slice is created as follows:
       Slice (nv, xyzv, [, val [, plane [, iso [, opt_3d]]]])
   where :
     nv is a one-dimensional integer array whose ith entry is the
        number of vertices of the ith face.
     xyzv is a two-dimensional array dimensioned sum (nv) by 3.
        The first nv [0] triples are the vertices of face [0],
        the next nv [1] triples are the vertices of face [1], etc.
     val (if present) is an array the same length as nv whose
        ith entry specifies a color for face i.
     plane (if present) says that this is a plane Slice, and
        all the vertices xyzv lie in this plane.
     iso (if present) says that this is the isosurface for
        the given value.
     opt_3d (if present) specifies the 3d options used to
        color the slice.
     mask (if present) specifies how the surface is to be masked.
     contours (if present) is either an array of contour values
        or an integer giving the number of contours.
     scale (if present) is the contour scale: "lin", "log",
        or "normal".
     edges (if present) is 1 for showing edges, 0 for not.
   A Slice object or two Slice objects are created by a call
   to the function sslice (q. v.). The function sslice accepts
   either a mesh and a specification of how to slice it
   (isosurface or plane), or else a previously created Slice,
   a plane to slice it with, and whether you want to keep
   the resulting "front" Slice or both Slices.
   """

   def type (self) :
       return Slice3dType

   def __init__ (self, nv, xyzv, val = None, plane = None, iso = None,
       smooth = None, opt_3d = ["wm", "f3"], mask = "max", contours = None,
       scale = "lin", edges = 0) :
       self.nv = nv
       self.xyzv = xyzv
       self.val = val
       self.plane = plane
       self.iso = iso
       self.smooth = smooth
       if type (opt_3d) != ListType :
          self.opt_3d = [opt_3d]
       else :
          self.opt_3d = opt_3d
       self.mask = mask
       self.contours = contours
       self.scale = scale
       self.edges = edges
       self.mesh_type = "xy"
       self.z_c_switch = 0
       self.z_contours_scale = scale
       self.c_contours_scale = scale
       if type (self.contours) == IntType or self.contours is None :
          self.number_of_z_contours = self.contours
          self.number_of_c_contours = self.contours
          self.z_contours_array = None
          self.c_contours_array = None
       elif type (self.contours) == ArrayType :
          self.z_contours_array = self.contours
          self.c_contours_array = self.contours
          self.number_of_z_contours = len (self.contours)
          self.number_of_c_contours = len (self.contours)
       if edges == 0 :
          self.edges = "w3" in self.opt_3d or \
                       "w4" in self.opt_3d or \
                       "wm" in self.opt_3d

   def __del__ (self) :
       del self.nv
       del self.xyzv
       if self.val is not None: del self.val
       if self.plane is not None: del self.plane
       del self.mask
       if self.iso is not None: del self.iso
       if self.smooth is not None: del self.smooth
       del self.opt_3d
       if self.contours is not None: del self.contours
       del self.scale
       del self.edges

   def new (self, nv, xyzv, val = None, plane = None, iso = None,
      smooth = None, opt_3d = ["wm", "f3"], mask = "max", contours = None,
      scale = "lin", edges = 0) :
       del self.nv
       del self.xyzv
       del self.val
       del self.plane
       del self.mask
       del self.iso
       del self.smooth
       del self.opt_3d
       del self.contours
       del self.scale
       del self.edges
       self.__init__ (nv, xyzv, val, plane, iso, smooth, opt_3d,
          mask, contours, scale, edges)

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       Slice characteristics. No error checking is done.
       """
       for k in keywords.keys ():
          setattr (self, k, keywords [k])
       if type (self.opt_3d) != ListType :
          self.opt_3d = [self.opt_3d]
       if self.edges == 0 :
          self.edges = "w3" in self.opt_3d or \
                       "w4" in self.opt_3d or \
                       "wm" in self.opt_3d

_SliceError = "SliceError"

def sslice ( v1, v2, varno = 1, nslices = 1, opt_3d = None) :
    """This function slices a Mesh3d or a previous Slice object,
    and returns a new Slice object, depending on how it is called.
       sslice (m, plane, varno = 1)
    returns a plane Slice through the specified mesh m. varno
    is the number of the variable to be used to color the Slice
    (1-based).
       sslice (m, val, varno = 1)
    returns an isosurface Slice of m. varno is the number of
    the variable with respect to which the Slice is taken (1-based)
    and val is its value on the isosurface. varno can be let default
    to 1 if there is only one function to worry about.
       sslice (s, plane, nslices = 1)
    slices the isosurface s and returns the part in "front" of the
    plane if nslices == 1, or returns a pair of Slices [front, back]
    if nslices == 2. If s is a Surface, then this operation will
    produce a Slice or pair of Slices from the given Surface.
    N. B. If you want just the "back" surface, you can achieve this
    by calling sslice with - plane instead of plane.
       The opt_3d argument can be used to overrule defaults, which
    are intelligently figured out, but may not always be what you want.
    """
    if v1.type () == Mesh3dType :
       if type (v1.c) != ListType :
          funcs = [v1.c]
       else :
          funcs = v1.c
       if not hasattr (v2, "type") or v2.type () != PlaneType :
          try :
             val = float (v2)
          except :
             raise _SliceError, \
                "second argument is not coercible to float."
          plane = None
       else :
          plane = v2
          val = None
       if opt_3d is None :
          if plane is not None:
             opt_3d = ["wm", "f4"]
          else :
             opt_3d = ["f4"]

       # First we build a "mesh3" if v1 does not already have one.
       if v1.structured :           #and not hasattr (v1, "m3") :
          setattr (v1, "m3", mesh3 (v1.x, v1.y, v1.z, funcs = funcs))
       elif not v1.structured and v1.avs == 1 : #and not hasattr (v1, "m3") :
          if no_of_dims (v1.x) == 3 and \
             no_of_dims (v1.y) == 3 and \
             no_of_dims (v1.z) == 3 :
             setattr (v1, "m3", mesh3 (array ( [v1.x, v1.y, v1.z], Float),
                funcs = funcs))
          elif no_of_dims (v1.x) == 1 and \
             no_of_dims (v1.y) == 1 and \
             no_of_dims (v1.z) == 1 :
             setattr (v1, "m3", mesh3 (v1.x, v1.y, v1.z, funcs = funcs,
                verts = v1.get_verts_list ()))
          else :
             raise _SliceError, \
                "x, y, and z must have either 1 or 3 dimensions."
       elif not hasattr (v1, "m3") :
          raise _SliceError, \
             "Option not implemented."
       # A Yorick 'mesh3' object now exists, slice it as specified.
       if plane is None : # An isosurface Slice
          [nv, xyzv, dum] = slice3 (v1.m3, varno, None, None, value = val)
          return Slice (nv, xyzv, iso = val, opt_3d = opt_3d)
       else :             # a plane Slice
          [nv, xyzv, val] = slice3 (v1.m3, plane.rep (), None, None, varno )
          return Slice (nv, xyzv, val, plane = plane, opt_3d = opt_3d)
    elif v1.type () == SurfaceType :
       plane = v2
       if not hasattr (v2, "type") or v2.type () != PlaneType :
          raise _SliceError, \
             "Second argument is not a Plane."

       if v1.z_c_switch :
          z = v1.c
          c = v1.z
          if "s3" in v1.opt_3d :
             scale = v1.c_contours_scale
             if v1.c_contours_array is not None:
                contours = v1.c_contours_array
             else :
                contours = v1.number_of_c_contours
          elif "s4" in v1.opt_3d :
             scale = v1.z_contours_scale
             if v1.z_contours_array is not None:
                contours = v1.z_contours_array
             else :
                contours = v1.number_of_z_contours
       else :
          z = v1.z
          c = v1.c
          if "s3" in v1.opt_3d :
             scale = v1.z_contours_scale
             if v1.z_contours_array is not None:
                contours = v1.z_contours_array
             else :
                contours = v1.number_of_z_contours
          elif "s4" in v1.opt_3d :
             scale = v1.c_contours_scale
             if v1.c_contours_array is not None:
                contours = v1.c_contours_array
             else :
                contours = v1.number_of_c_contours
       edges = "w3" in v1.opt_3d or "w4" in v1.opt_3d or "wm" in v1.opt_3d \
            or "i3" in v1.opt_3d or "i4" in v1.opt_3d
       if "s4" in v1.opt_3d :
          [nv, xyzv, col] = slice3mesh (v1.x, v1.y, z, color = c, smooth = 1)
       elif "f3" in v1.opt_3d :
          [nv, xyzv, col] = slice3mesh (v1.x, v1.y, z, color = z, smooth = 0)
          contours = None
          scale = "lin"
       elif "f4" in v1.opt_3d :
          [nv, xyzv, col] = slice3mesh (v1.x, v1.y, z, color = c, smooth = 0)
          contours = None
          scale = "lin"
       else :
          [nv, xyzv, col] = slice3mesh (v1.x, v1.y, z)
          contours = None
          scale = "lin"
       if nslices == 1 :                                    
          [nv, xyzv, col] = slice2 (plane.rep (), nv, xyzv, col)
          return Slice (nv, xyzv, col, opt_3d = v1.opt_3d, mask = v1.mask,
             contours = contours, scale = scale, edges = edges)
       elif nslices == 2 :
          [nv, xyzv, col, nb, xyzb, valb] = slice2x (plane.rep (),
              nv, xyzv, col)
          return [Slice (nv, xyzv, col, opt_3d = v1.opt_3d, mask = v1.mask,
             contours = contours, scale = scale, edges = edges),
                  Slice (nb, xyzb, valb, opt_3d = v1.opt_3d, mask = v1.mask,
             contours = contours, scale = scale, edges = edges)]
       else :
          raise _SliceError, "Illegal number (" + `nslices` + \
             ") of slices requested."
    elif v1.type () == Slice3dType :
       plane = v2
       if not hasattr (v2, "type") or v2.type () != PlaneType :
          raise _SliceError, \
             "Second argument is not a Plane."

       if nslices == 1 :
          [nv, xyzv, val] = slice2 (plane.rep (), v1.nv, v1.xyzv,
              v1.val)
          return Slice (nv, xyzv, val, plane = v1.plane, iso = v1.iso,
                        opt_3d = v1.opt_3d)
       elif nslices == 2 :
          [nf, xyzf, valf, nb, xyzb, valb] = slice2x (plane.rep (), v1.nv,
              v1.xyzv, v1.val)
          return [Slice (nf, xyzf, valf, plane = v1.plane, iso = v1.iso,
                         opt_3d = v1.opt_3d),
                  Slice (nb, xyzb, valb, plane = v1.plane, iso = v1.iso,
                         opt_3d = v1.opt_3d)]
       else :
          raise _SliceError, "Illegal number (" + `nslices` + \
             ") of slices requested."
    else :
       raise _SliceError, "Can only slice a Mesh3d or another Slice."
