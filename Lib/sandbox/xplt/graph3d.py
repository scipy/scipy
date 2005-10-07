# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# Nar, Numeric, and shapetest are all imported by graph:

from types import *
from graph import *
from graftypes import *

class Graph3d ( Graph ) :

   """
   g = Graph3d ( <surface  list>, ...keyword arguments...) will create
   a three-dimensional graphics object consisting of several surfaces
   plus a global environment for the object. It will accept one or
   a list of Plotter objects or plotter_identifiers, or will try
   to complete a generic connection of its own if asked to plot
   without such a plotter specification.
   <surface list> is one or a sequence of Surface  objects.
   The keyword arguments for Graph3d are:

       plotter = <Plotter object> or a sequence of <Plotter object>s
              if you have a plotter object or a sequence of them that
              you want the curve to use when plotting itself. I
              recommend against this; a curve can create its own
              plotter if you don't give it one.
       filename = <string> or a sequence of <string>s if you want to
              connect to Narcisse with the named file(s). (The default
              is " ".) Only one of the keywords 'plotter' or 'filename'
              is allowed, and both are optional.
              NOTE: the possibility of a sequence of file names or
              plotters allows one to display Narcisse graphs on one's
              own machine as well as one or more remote machines.
              In Gist, the filename argument is the same as the
              display argument (below).
       display = <string> or a sequence of <string>s if you want to
              display on the named hosts. The form of <string> is
              the usual "hostname:server.screen". The purpose of
              this argument is to allow you to continue without
              exiting python if you have to open a Gist window
              without the DISPLAY environment variable having been
              set, or if you want to open Gist windows on more
              than one host.
       titles = <value> where <value> is a string or a sequence
                of up to four strings, giving the titles in the
                order bottom, top, left, right.
       title_colors = <value> where value is an integer or string
                or a sequence of up to four integers or strings giving
                the colors of the titles.
       grid_type = <string> where "none" means no axis grid;
                   "axes" means a pair of axes with tick marks;
                   "wide" means a widely spaced 2d grid; and
                   "full" means a closely spaced 2d grid.
       axis_labels = <value> where <value> is a string or
                   sequence of up to five strings representing
                   the labels of the x axis, the y axis, the z
                   axis, the c axis, and the right y axis.
                   For Gist, only the first three count, and
                   if necessary, will be truncated to a single
                   character.
       gnomon = (0|1) (Gist only). If 1, show the gnomon (a simple display
                   showing the orientation of the xyz axes of the object).
       x_axis_label, y_axis_label, z_axis_label, c_axis_label, and
                   yr_axis_label may be used to change individual
                   axis labels.
       axis_limits = <value> where <value> is a pair [xmin, xmax] or
                   a sequence of up to five pairs, where the second
                   would be the y limits, the third the z limits,
                   the fourth the c limits, and the fifth the right
                   y limits.
       x_axis_limits, y_axis_limits, z_axis_limits, c_axis_limits, and
                   yr_axis_limits may be used to change individual
                   axis limits.
       x_factor, y_factor (defaults 1.0). A factor > 1.0 compresses
                   the graph by that factor in the specified direction
                   (relative to the plane of the graph). A factor < 1.0
                   would enlarge it in the specified direction. So far
                   available only in Gist.
       axis_scales = <value> where <value> is a string or a sequence
                   of up to five strings, each of which is either "lin"
                   (linear scale) or "log" (logarithmic scale). Omitted
                   strings default to "lin". The order thay are
                   specified is x, y, z, c, and right y. If this
                   argument is missing, the existing axis scales
                   are unchanged.
       x_axis_scale, y_axis_scale, z_axis_scale, c_axis_scale, and
                   yr_axis_scale may be used to change individual
                   axis scales.
       text = <value> where value is one or a sequence of strings
                   representing texts to be placed on the plot.
       text_color = <value> where <value> is one or a sequence
                   of integer color numbers or strings (names of
                   common colors) giving colors for the texts.
       text_size = <value> where <value> is one or a sequence of
                   integers giving (roughly) the number of characters
                   in a line on the graph. The larger this number, the
                   smaller the size.
       text_pos = <value> where <value> is a pair or a sequence of
                   reals between 0. and 1.0 giving the relative
                   position of the lower left corner of a text
                   in the graphics window.
       Narcisse viewing angles:
          phi = <integer value> specifies the angle that the line from
                      the view point to the origin makes with the positive
                      z axis. The angle is in degrees. (default 45)
          theta = <integer value> specifies the angle made by the projection
                      of the line of view on the xy plane with the
                      positive x axis. The angle is in degrees. (default 45)
          roll = <integer value> specifies the angle of rotation of the
                      graph around the line from the origin to the view point.
                      The angle is measured in the positive direction,
                      i. e., if your right thumb is aligned outwards along
                      the line from the origin to the view point, then
                      your fingers curl in the positive direction.
                      The angle is in degrees. Not implemented for
                      Gist graphics. (default 0)
       Gist viewing angles:
          theta = <integer value> specifies the angle that the z axis
                makes with its projection onto the surface of the graph;
                positive if the z axis points out of the screen,
                negative if it points inwards. (default 30)
                If the z axis is in the screen, then the x and y axes
                will appear horizontal. Note: If the gnomon is turned on,
                then axes in the plane of the graphics window or pointing
                into it will have their labels highlighted.
          phi = <integer value> specifies the angle through which the
                x axis is rotated in the xy plane; positive if in the
                direction of the y axis, negative otherwise. (default -45)
       distance = <integer value> specifies the distance of the view
                   point from the origin. This is an integer between
                   0 and 20. 0 makes the distance infinite; otherwise
                   the smaller the number, the closer you are. This
                   number does not affect the size of the graph, but
                   rather the amount of distortion (on account of
                   perspective) in the picture (the closer you are,
                   the more distortion).
       link = <value> Used to link surfaces of different 3d options.
           normally all surfaces in a graph will have the same
           3d options. This value should be set to 1 if you want to
           graph two or more surfaces with different 3d options.
           otherwise multiple surface graphs will appear with the
           options of the last surface specified.
           Narcisse only.
       connect = <value> set to 1 for graphs of more than one surface
           to provide better hidden line removal. Must not be used
           with link.
           Narcisse only.
       sync = <value> set to 1 to synchronize with Narcisse before
           plotting the next graph. Keeps graphs sent in rapid
           succession from becoming garbled. Defaults to 1; set to
           0 if you don't have a timing problem.
           Narcisse only.
       Lighting variables (so far available only in Gist):
       ambient = <value> is a light level (arbitrary units)
           that is added to every surface independent of its orientation.
       diffuse = <value> is a light level which is proportional
           to cos(theta), where theta is the angle between the surface
           normal and the viewing direction, so that surfaces directly
           facing the viewer are bright, while surfaces viewed edge on are
           unlit (and surfaces facing away, if drawn, are shaded as if they
           faced the viewer).
       specular = <value> 
       spower = <value>
       sdir = <value> 
           specular = S_LEVEL is a light level proportional to a high
           power spower=N of 1+cos(alpha), where alpha is the angle between
           the specular reflection angle and the viewing direction.  The light
           source for the calculation of alpha lies in the direction XYZ (a
           3 element vector) in the viewer's coordinate system at infinite
           distance.  You can have ns light sources by making S_LEVEL, N, and
           XYZ (or any combination) be vectors of length ns (3-by-ns in the
           case of XYZ).
       z_scale = <real value>
           if not equal to 1, affects the depth in z.
       color_bar = 0 or 1 (1 enables plotting of a color bar on
                any graphs for which it is meaningful (colored contour
                plots, filled contour plots, cell arrays, filled
                meshes and polygons).
       color_bar_pos (ignored unless a color bar is actually plotted)
                is a 2d array [ [xmin, ymin], [xmax, ymax]] specifying
                where (in window coordinates) the diagonally opposite
                corners of the color bar are to be placed.
                (the z_contours_array and/or c_contours array can be used
                to control precisely what is shown in the color bar.)
                currently not implemented in Gist.
       split = 0 or 1 (default 1) If on, causes the palette to be split
                when both planes and isosurfaces are present in a graph,
                so that isosurfaces are shaded according to current
                light settings, while plane sections of the mesh are
                colored according to a specified function.
                Currently not implemented in Narcisse.
   """

   def type (self) :
       return Graph3dType

   _color_card_dict = { "absolute" : 0 , "binary" : 1 ,
   "bluegreen" : 2 , "default" : 3 , "negative" : 4 , "positive" : 5 ,
   "rainbow" : 6 , "rainbowhls" : 7 , "random" : 8 , "redblue" : 9 ,
   "redgreen" : 10 , "shifted" : 11 }

   _axes = ["x", "y", "z", "c", "yr"]

   _NotASurface = "NotASurface"
   _LinkError = "LinkError"
   _SurfaceSpecError = "SurfaceSpecError"
   _SurfaceChangeError = "SurfaceChangeError"

   # The following are actually used in Graph
   _limits_keywords = ["x_axis_limits", "y_axis_limits", "z_axis_limits",
                       "c_axis_limits", "yr_axis_limits"]
   _scale_keywords = ["x_axis_scale", "y_axis_scale", "z_axis_scale",
                      "c_axis_scale", "yr_axis_scale"]
   _label_keywords = ["x_axis_label", "y_axis_label", "z_axis_label",
                      "c_axis_label", "yr_axis_label"]
   _no_of_axes = 5

   def __init__ ( self , surface_list , *kwds , ** keywords ) :
      if len ( kwds ) == 1 :
         keywords = kwds[0]
      if is_scalar ( surface_list ) :
         if surface_list.type () != SurfaceType and \
            surface_list.type () != Mesh3dType and \
            surface_list.type () != Slice3dType :
           raise self._NotASurface , \
           "Attempt to initialize with something that is not a surface, " + \
            "slice,  or mesh."
         self._s = [surface_list]
      else :
         for i in range ( len ( surface_list ) ) :
            if surface_list[i].type () != SurfaceType and \
               surface_list[i].type () != Mesh3dType and \
               surface_list[i].type () != Slice3dType :
               raise self._NotASurface , \
               "Attempt to initialize with something that is not a surface, " \
               + "slice,  or mesh."
         self._s = surface_list
      self._s_ln = len ( self._s )
      if keywords.has_key ("x_factor") :
         self._x_factor = keywords ["x_factor"]
      else :
         self._x_factor = 1.0
      if keywords.has_key ("y_factor") :
         self._y_factor = keywords ["y_factor"]
      else :
         self._y_factor = 1.0
      if keywords.has_key ("link") :
         if keywords.has_key ("connect") :
            raise _LinkError , "Can't specify both link and connect."
         self._link = keywords ["link"]
      else :
         self._link = 0
      if keywords.has_key ("connect") :
         self._connect = keywords ["connect"]
      else :
         self._connect = 0
      if keywords.has_key ("phi") :
         self._phi = keywords ["phi"]
      else :
         self._phi = None
      if keywords.has_key ("theta") :
         self._theta = keywords ["theta"]
      else :
         self._theta = None
      if keywords.has_key ("roll") :
         self._roll = keywords ["roll"]
      else :
         self._roll = None
      if keywords.has_key ("distance") :
         self._distance = keywords ["distance"]
      else :
         self._distance = 0
      if keywords.has_key ("z_scale") :
         self.z_scale = keywords ["z_scale"]
      else :
         self.z_scale = 1.0
      if keywords.has_key ("gnomon") :
         self._gnomon = keywords ["gnomon"]
      else :
         self._gnomon = None
      if keywords.has_key ( "color_card" ) :
         self._color_card = keywords ["color_card"]
      else :
         self._color_card = "default"
      self._axis_limits = [[0., 0.], [0., 0.], [0., 0.], [0., 0.], [0., 0.]]
      self._axis_scales = ["lin", "lin", "lin", "lin", "lin"]
      self._axis_labels = ["X axis", "Y axis", "Z axis", "C axis", "YR axis"]
      # Gist lighting
      self.lighting_dict = {}
      if keywords.has_key ("ambient") :
         self.lighting_dict ["ambient"] = keywords ["ambient"]
      else :
         self.lighting_dict ["ambient"] = None
      if keywords.has_key ("diffuse") :
         self.lighting_dict ["diffuse"] = keywords ["diffuse"]
      else :
         self.lighting_dict ["diffuse"] = None
      if keywords.has_key ("specular") :
         self.lighting_dict ["specular"] = keywords ["specular"]
      else :
         self.lighting_dict ["specular"] = None
      if keywords.has_key ("spower") :
         self.lighting_dict ["spower"] = keywords ["spower"]
      else :
         self.lighting_dict ["spower"] = None
      if keywords.has_key ("sdir") :
         self.lighting_dict ["sdir"] = keywords ["sdir"]
      else :
         self.lighting_dict ["sdir"] = None
      if keywords.has_key ("style") :
         self._style = keywords ["style"]
      else :
         self._style = "nobox.gs"
      if keywords.has_key ("grid_type") :
         self._grid_type = keywords ["grid_type"]
      else :
         self._grid_type = "none"
      if keywords.has_key ( "opt_3d" ) :
         self.opt_3d = keywords ["opt_3d"]
      else :
         self.opt_3d = "wm"
      if keywords.has_key ( "mask" ) :
         self.mask = keywords ["mask"]
      else :
         self.mask = "max"
      if keywords.has_key ( "split" ) :
         self._split = keywords ["split"]
      else :
         self._split = 0
      # Everything else is Graph generic:
      self._label_type = " "
      Graph.__init__ ( self , keywords )

   def new ( self , surface_list , ** keywords ) :
      """new ( surface_list, <keyword arguments> ) cleans out a Graph3d
      and reinitializes it. This has the same argument list as
      Graph3d. Do not change the plotter list or filename list (to
      avoid the tedious on-and-off flickering of windows) unless
      this is actually requested.
      """
      pl = self._plotter_list
      fl = self._filename_list
      self._plotter_list = []
      self._filename_list = []
      del self._s, self._s_ln, self._link, self._connect, self._phi, \
          self._theta, self._roll, self._distance, self._titles, \
          self._title_colors, self._text, self._text_color, self._text_size, \
          self._sync, self._grid_type, self._axis_labels, self._axis_scales, \
          self._axis_limits
      self.__init__ ( surface_list , keywords )
      if self._plotter_list == [] :
         self._plotter_list = pl
         self._filename_list = fl

   def add ( self , surface ) :
      """ add ( surface ) adds a surface with its characteristics to the
      existing plot. Surfaces are numbered in the order that they are
      added, beginning with 1.
      """
      if surface.type () != SurfaceType and surface.type () != Mesh3dType and \
         surface_list[i].type () != Slice3dType :
         raise self._NotASurface , \
         "Attempt to add something not a surface, slice, or mesh."
      self._s.append ( surface )
      self._s_ln = len (self._s)

   def delete ( self , n ) :
      """delete ( n ) deletes the nth surface from the Graph. Note
      that surfaces are numbered beginning with 1.
      """
      DeleteError = "DeleteError"
      if 1 <= n <= self._s_ln :
         self._s [n-1:n] = []
         self._s_ln = len (self._s)
      else :
         raise DeleteError , "There is no surface numbered " + `n` + \
            " in the current graph, which has only " + `self._s_ln` + \
            " surfaces."

   def set_surface_list (self, surface_list) :
      del self._s
      self._s = surface_list
      self._s_ln = len (surface_list)
      return

   def replace ( self , n , surface ) :
      """replace ( n , surface ) replaces the nth surface in the Graph
      with the specified surface. Note that surfaces are numbered
      beginning with 1.
      """
      if surface.type () != SurfaceType and surface.type () != Mesh3dType and \
               surface.type () != Slice3dType :
         raise self._NotASurface , \
            "Attempt to add something not a surface, slice, or mesh."
      ReplaceError = "ReplaceError"
      if 1 <= n <= self._s_ln :
         self._s [n-1] = surface
      else :
         raise ReplaceError , "There is no surface numbered " + `n` + \
            " in the surfacent graph, which has only " + `self._s_ln` + \
            " surfaces."

   def change_plot ( self , ** keywords ) :
      """change_plot ( <keyword arguments> ) is used to change any Graph3d
      characteristics except the surfaces being graphed. Use the add,
      delete, and/or replace commands to do that. change_plot will
      draw the graph without sending surface coordinates, unless
      keyword send is 1. Generally, change_plot should be used when
      the graph needs to be recomputed, and quick_plot when it does not.
      change_plot does no error checking and does not convert user-friendly
      names of colors and such into numbers.
      """
      for k in keywords.keys ():
         if k == "surface" or k == "mesh" :
            raise self._SurfaceChangeError, \
               "Use add, delete, or replace to change surfaces in a graph."
         setattr (self, "_" + k, keywords [k])
      if "send" in keywords.keys ():
         send = keywords ["send"]
      else:
         send = 0
      self._send_coordinates = send
      self.plot ( )
      self._send_coordinates = 1

   def quick_plot ( self , ** keywords ) :
      """quick_plot ( <keyword arguments> ) is used to change some Graph3d
      characteristics which do not demand that the graph be recomputed.
      You can change the characteristics of a surface in the graph by
      specifying its number (surface = n) and any combination of the
      traits color_card, opt_3d, mesh_type, or mask. Or you can change
      such overall graph characteristics as  titles, title_colors, text,
      text_color, text_size, text_pos, color_card, grid_type, sync, theta,
      phi, roll, and axis_labels.
      Or you can do both.
      The changes will be effected and the graph redrawn.
      Things that you cannot change include axis limits and scales,
      and the coordinates of a surface. Use change_plot if axis limits
      and scales are among the things you want to change, and use add,
      delete, or replace followed by a call to plot, if you wish to
      change a surface.
      quick_plot will not work right for link'ed surfaces. Once the
      changes have been made, you will have to call plot.
      """
      ChangeError = "ChangeError"
      PlottersNotStarted = "PlottersNotStarted"
      if not self._plotters_started :
         raise PlottersNotStarted , \
            "quick_plot requires that all plotters have already been started."
      if keywords.has_key ("opt_3d") :
         self._opt_3d = keywords ["opt_3d"]
         self.opt_3d_change = 1
      if keywords.has_key ("mesh_type" ) :
         self._mesh_type = keywords ["mesh_type"]
         self.mesh_type_change = 1
      if keywords.has_key ("mask") :
         self._mask = keywords ["mask"]
         self.mask_change = 1
      if keywords.has_key ("color_card") :
         if type ( keywords [ "color_card" ] ) == IntType :
            self._color_card = keywords [ "color_card" ]
            self.color_card_change = 1
         elif self._color_card_dict.has_key \
            ( keywords [ "color_card" ] ) :
            self._color_card = \
               self._color_card_dict [keywords [ "color_card" ]]
            self.color_card_change = 1
      if keywords.has_key ( "surface" ) or keywords.has_key ( "mesh" ) :
         if keywords.has_key ( "surface" ) :
            self.n = keywords ["surface"]
         else :
            self.n = keywords ["mesh"]
         self.color_card_change = 0
         self.opt_3d_change = 0
         self.mesh_type_change = 0
         self.mask_change = 0
         if 1 <= self.n <= self._s_ln :
            if keywords.has_key ("color_card") :
               if type ( keywords [ "color_card" ] ) == IntType :
                  self._s [self.n - 1].color_card = keywords [ "color_card" ]
                  self.color_card_change = 1
               elif self._color_card_dict.has_key \
                  ( keywords [ "color_card" ] ) :
                  self._s [self.n - 1].color_card = \
                     self._color_card_dict [keywords [ "color_card" ]]
                  self.color_card_change = 1
         else :
            raise ChangeError , "There is no surface numbered " + `n` + \
            " in the current graph, which has only " + `self._s_ln` + \
            " surfaces."
      if keywords.has_key ("grid_type") :
         self._grid_type = keywords ["grid_type"]
      if keywords.has_key ("sync") :
         self._sync = keywords ["sync"]
      if keywords.has_key ("theta") :
         self._theta = keywords ["theta"]
      if keywords.has_key ("phi") :
         self._phi = keywords ["phi"]
      if keywords.has_key ("roll") :
         self._roll = keywords ["roll"]
      if keywords.has_key ("ambient") :
         self.lighting_dict ["ambient"] = keywords ["ambient"]
      else :
         self.lighting_dict ["ambient"] = None
      if keywords.has_key ("diffuse") :
         self.lighting_dict ["diffuse"] = keywords ["diffuse"]
      else :
         self.lighting_dict ["diffuse"] = None
      if keywords.has_key ("specular") :
         self.lighting_dict ["specular"] = keywords ["specular"]
      else :
         self.lighting_dict ["specular"] = None
      if keywords.has_key ("spower") :
         self.lighting_dict ["spower"] = keywords ["spower"]
      else :
         self.lighting_dict ["spower"] = None
      if keywords.has_key ("sdir") :
         self.lighting_dict ["sdir"] = keywords ["sdir"]
      else :
         self.lighting_dict ["sdir"] = None
      Graph.change ( self , keywords ) # Change generic traits
      for ipl in range (len (self._plotter_list)) :
         pl = self._plotter_list [ipl]
         pl.quick_plot (self)

   def plot ( self ) :
       """plot ( ) plots a 3d graph object. If the user has not by
       now specified plotter(s) or filename(s) then a generic plotter
       object will be created, if it is possible to find a local
       Graphics routine.
       """
       self._init_plotter ( )
       # grid_type means something different on 3d graphs
       for ipl in range (len (self._plotter_list)) :
          pl = self._plotter_list [ipl]
          pl.plot3d (self)
   
   def move_light_source (self, ** keywords) :
       """move_light_source (nframes = 30, angle = 12) causes a
       movie of the surrent graph to be shown, with the light source
       rotating through angle degrees each frame, for a total
       of nframes. If angle is not specified and nframes is,
       then the angle defaults to 360 / nframes.
       """
       if not keywords.has_key ("angle") and keywords.has_key ("nframes") :
          nframes = keywords ["nframes"]
          angle = 2 * pi / nframes
       else :
          if keywords.has_key ("nframes") :
             nframes = keywords ["nframes"]
          else :
             nframes = 30
          if keywords.has_key ("angle") :
             angle = keywords ["angle"]
          else :
             angle = 2 * pi / nframes
       for ipl in range (len (self._plotter_list)) :
          pl = self._plotter_list [ipl]
          pl.move_light_source (self, angle, nframes)

   def rotate (self, ** keywords) :
       """rotate (axis = array ([-1., 1., 0.],  Float), angle = 12,
       nframes = 30) Causes the current graph to be rotated about
       the specified axis, through the specified angle each frame,
       for a total of nframes frames. If angle is not specified and
       nframes is, then the angle defaults to 360 / nframes.
       """
       if keywords.has_key ("axis") :
          axis = keywords ["axis"]
       else :
          axis = array ([-1., 1., 0.],  Float)
       if not keywords.has_key ("angle") and keywords.has_key ("nframes") :
          nframes = keywords ["nframes"]
          angle = 2 * pi / nframes
       else :
          if keywords.has_key ("nframes") :
             nframes = keywords ["nframes"]
          else :
             nframes = 30
          if keywords.has_key ("angle") :
             angle = pi * keywords ["angle"] / 180. # angle is degrees
          else :
             angle = 2 * pi / nframes
       for ipl in range (len (self._plotter_list)) :
          pl = self._plotter_list [ipl]
          pl.rotate_graph (self, axis, angle, nframes)
