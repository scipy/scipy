# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from graph3d import *
from mesh3d import *
from surface import *

from Numeric import *
from scipy_base.fastumath import *

from arrayfns import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

def demo (*i) :
   x = span (-1, 1, 64, 64)
   y = transpose (x)
   z = (x + y) * exp (-6.*(x*x+y*y))

   if len (i) == 0 or i [0] == 1 :
      s1 = Surface (z = z, opt_3d = "wm", mask = "sort")

      g1 = Graph3d (s1, color_card = "gray.gp", titles = "opaque wire mesh", y_factor = 2.)

      g1.plot ()

      paws ()

      s1.set (mask = "none")

      g1.change (titles = "transparent wire mesh")
      g1.plot ()
      paws ()

      s1.set (ecolor = "red")
      g1.change (titles = "transparent wire mesh in red")
      g1.plot ()
      paws ()

      s1.set (mask = "sort", shade = 1)
      g1.change (titles = "opaque shaded mesh with red lines")
      g1.plot ()
      paws ()

      s1.set (opt_3d = "none")
      g1.change (titles = "opaque shaded mesh with no lines")
      g1.plot ()
      paws ()

      g1.change (titles = "same with different lighting")
      g1.quick_plot (diffuse=.1, specular = 1., sdir=array([0,0,-1]))
      paws ()
   
      s1.set (c = z - min (ravel (z)), opt_3d = "s4", shade = 0)
      g1.change (color_card = "rainbow.gp",
         titles = "surface shaded according to height")
      g1.plot ()
      paws ()

   if len (i) == 0 or i [0] == 2 :
      s1 = Surface (z = z, opt_3d = "none", mask = "sort", shade = 1)

      # The purpose of the try clause is to prevent creation of a new Graph
      # object if g1 already exists. Creation of a new one forces the
      # current Plotter(s) to close and new one(s) to open. When working
      # with Xwindows, this leads to annoying popping of windows in and
      # out of existence.
      try:
         g1 = g1
      except:
         g1 = Graph3d (s1, ambient = 0.2, diffuse = .2, specular = 1.,
            color_card = "gray.gp", titles = "backlit surface", y_factor = 2.)
      else:
         g1.new (s1, ambient = 0.2, diffuse = .2, specular = 1.,
            color_card = "gray.gp", titles = "backlit surface", y_factor = 2.)

      g1.plot ()
      paws ()
      g1.change (titles = "moving light source")
      g1.move_light_source ()
      paws ()

   if len (i) == 0 or i [0] == 3 :
      nx = 20
      ny = 20
      nz = 20
      xyz = zeros ( (3, nx, ny, nz), Float)
      xyz [0] = multiply.outer ( span (-1, 1, nx), ones ( (ny, nz), Float))
      xyz [1] = multiply.outer ( ones (nx, Float),
         multiply.outer ( span (-1, 1, ny), ones (nz, Float)))
      xyz [2] = multiply.outer ( ones ( (nx, ny), Float), span (-1, 1, nz))
      r = sqrt (xyz [0] ** 2 + xyz [1] **2 + xyz [2] **2)
      theta = arccos (xyz [2] / r)
      phi = arctan2 (xyz [1] , xyz [0] + logical_not (r))
      y32 = sin (theta) ** 2 * cos (theta) * cos (2 * phi)

      m1 = Mesh3d (x = span (-1, 1, nx),
                   y = span (-1, 1, ny),
                   z = span (-1, 1, nz),
                   c = r * (1. + y32))
      del r, theta, phi, xyz, y32

      s1 = sslice (m1, .50, varno = 1) # , opt_3d = "f3") # (inner isosurface)
      s2 = sslice (m1, 1.0, varno = 1) # , opt_3d = "f3") # (outer isosurface)

      pxy = Plane (array ([0., 0., 1.], Float ), zeros (3, Float))
      pyz = Plane (array ([1., 0., 0.], Float ), zeros (3, Float))

      # pseudo-colored plane slice, then cut it in half and save
      # only the front half
      s3 = sslice (m1, pyz, opt_3d = ["wm", "f4"])
      s3 = sslice (s3, pxy, nslices = 1, opt_3d = ["wm", "f4"])

      # cut the inner isosurface in half so that we can slice the
      # top off one of the halves:
      [s1, s4] = sslice (s1, pxy, nslices = 2) #, opt_3d = "f3")
      s1 = sslice (s1, - pyz) #, opt_3d = "f3")

      # do the same with the outer isosurface:
      [s2, s5] = sslice (s2, pxy, nslices = 2) #, opt_3d = "f3")
      s2 = sslice (s2, - pyz) #, opt_3d = "f3")
     
      # The purpose of the try clause is to prevent creation of a new Graph
      # object if g1 already exists. Creation of a new one forces the
      # current Plotter(s) to close and new one(s) to open. When working
      # with Xwindows, this leads to annoying popping of windows in and
      # out of existence.
      try:
         g1 = g1
      except:
         g1 = Graph3d ([s3, s1, s4, s2, s5], gnomon = 1, color_card =
            "rainbowhls", diffuse = .2, specular = 1,
            mask = "min", split = 1) # opt_3d = ["f4"],
      else:
         g1.new ([s3, s1, s4, s2, s5], gnomon = 1, color_card =
            "rainbowhls", diffuse = .2, specular = 1,
            mask = "min", split = 1) 

      paws ()
      g1.plot ()
      paws ()
      g1.change (x_factor = 2., y_factor = 2.)
      g1.rotate ()
