# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from plwf import *
from pl3d import *
from movie import *
from slice3 import *
from yorick import *
from gist import *
from arrayfns import *
from MLab import rand

window3 (hcp = "talk.ps", dump = 1)

palette ("gray.gp")

demo5_n = 20. * ones (3)

making_movie = 0

def demo5_light (i) :
   global making_movie
   if i >= 30 : return 0
   theta = pi / 4 + (i - 1) * 2 * pi/29
   light3 (sdir = array ( [cos(theta), .25, sin(theta)], Float))
   # without an explicit call to draw3, the light3 function would
   # cause no changes until Python paused for input from the keyboard,
   # since unlike the primitive plotting functions (plg, plf, plfp, ...)
   # the fma call made by the movie function will not trigger the
   # 3D display list
   # any movie frame display function which uses the 3D drawing
   # functions in pl3d.i will need to do this
   # the !making_movie flag supresses the fma in draw3 if this function
   # is called by movie (which issues its own fma), but allows it
   # otherwise 

   draw3 ( not making_movie )
   return 1

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

from PR import *
def demo5 (*itest) :
   """demo5 () or demo5 (i)
     Run examples of use of pl3d.i, plwf.i, and slice3.i.  With
     argument I = 1, 2, or 3, run that particular demonstration.
     Read the source code to understand the details of how the
     various effects are obtained.

     demo5 (1) demonstrates the various effects which can be obtained
     with the plwf (plot wire frame) function.
     demo5 (2) demonstrates shading effects controlled by the light3
     function
     demo5 (3) demonstrates the slice3, slice2, and pl3tree functions,
     as well as changing the orientation of the 3D object
   """
   global making_movie
   if len (itest) == 0 or itest [0] == 1 :
      set_draw3_ (0)
      x = span (-1, 1, 64, 64)
      y = transpose (x)
      z = (x + y) * exp (-6.*(x*x+y*y))
      limits_(square = 1)
      print  "(plot wire frame) plwf,z,y,x"
      orient3 ( )
      light3 ( )
      plwf (z, y, x)
      [xmin, xmax, ymin, ymax] = draw3(1) # not necessary interactively
      limits (xmin, xmax, ymin, ymax)
      plt("opaque wire mesh", .30, .42)
      paws ( )
      print "plwf,z,y,x, shade=1,ecolor=\"red\""
      plwf(z,y,x,shade=1,ecolor="red")
      [xmin, xmax, ymin, ymax] = draw3(1) # not necessary interactively
      limits (xmin, xmax, ymin, ymax)
      paws()
      print "plwf,z,y,x, shade=1,edges=0"
      plwf(z,y,x,shade=1,edges=0)
      [xmin, xmax, ymin, ymax] = draw3(1) # not necessary interactively
      limits (xmin, xmax, ymin, ymax)
      paws ( )
      light3 ( diffuse=.1, specular=1., sdir=array([0,0,-1]))
      [xmin, xmax, ymin, ymax] = draw3(1)
      limits (xmin, xmax, ymin, ymax)
      paws ( )
      light3 ( diffuse=.5, specular=1., sdir=array([1,.5,1]))
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ( )
      light3 ( ambient=.1,diffuse=.1,specular=1.,
             sdir=array([[0,0,-1],[1,.5,1]]),spower=array([4,2]))
      [xmin, xmax, ymin, ymax] = draw3(1)
      limits (xmin, xmax, ymin, ymax)
      paws ( )
   if len (itest) == 0 or itest [0] == 2 :
      set_draw3_ (0)
      x = span (-1, 1, 64, 64)
      y = transpose (x)
      z = (x + y) * exp (-6.*(x*x+y*y))
      print "light3 function demo- default lighting"
      orient3 ( )
      light3 ( )
      plwf (z,y,x,shade=1,edges=0)
      [xmin, xmax, ymin, ymax] = draw3 (1) # not necessary interactively
      limits (xmin, xmax, ymin, ymax)
      paws( )
      print "light3,diffuse=.2,specular=1"
      light3(diffuse=.2,specular=1)
      limits_(square = 1)
      [xmin, xmax, ymin, ymax] = draw3(1) # not necessary interactively
      limits (xmin, xmax, ymin, ymax)
      paws()
      print "light3,sdir=[cos(theta),.25,sin(theta)]  -- movie"
      making_movie = 1
      movie(demo5_light, lims = [xmin, xmax, ymin, ymax])
      making_movie = 0
      fma()
      demo5_light(1)
      paws()
      light3()
   if len (itest) == 0 or itest [0] == 3 :
      nx = demo5_n [0]
      ny = demo5_n [1]
      nz = demo5_n [2]
      xyz = zeros ( (3, nx, ny, nz), Float)
      xyz [0] = multiply.outer ( span (-1, 1, nx), ones ( (ny, nz), Float))
      xyz [1] = multiply.outer ( ones (nx, Float),
         multiply.outer ( span (-1, 1, ny), ones (nz, Float)))
      xyz [2] = multiply.outer ( ones ( (nx, ny), Float), span (-1, 1, nz))
      r = sqrt (xyz [0] ** 2 + xyz [1] **2 + xyz [2] **2)
      theta = arccos (xyz [2] / r)
      phi = arctan2 (xyz [1] , xyz [0] + logical_not (r))
      y32 = sin (theta) ** 2 * cos (theta) * cos (2 * phi)
      m3 = mesh3 (xyz, funcs = [r * (1. + y32)])
      del r, theta, phi, xyz, y32

      print "   test uses " + `(nx - 1) * (ny - 1) * (nz - 1)` + " cells"
      elapsed = [0., 0., 0.]
      elapsed = timer_ (elapsed)
      elapsed0 = elapsed

      [nv, xyzv, dum] = slice3 (m3, 1, None, None, value = .50)
          # (inner isosurface)
      [nw, xyzw, dum] = slice3 (m3, 1, None, None, value = 1.)
          # (outer isosurface)
      pxy = plane3 ( array ([0, 0, 1], Float ), zeros (3, Float))
      pyz = plane3 ( array ([1, 0, 0], Float ), zeros (3, Float))
      [np, xyzp, vp] = slice3 (m3, pyz, None, None, 1)
          # (pseudo-colored slice)
      [np, xyzp, vp] = slice2 (pxy, np, xyzp, vp)
          # (cut slice in half)
      [nv, xyzv, d1, nvb, xyzvb, d2] = \
          slice2x (pxy, nv, xyzv, None)
      [nv, xyzv, d1] = \
          slice2 (- pyz, nv, xyzv, None)
          # (...halve one of those halves)
      [nw, xyzw, d1, nwb, xyzwb, d2] = \
          slice2x ( pxy , nw, xyzw, None)
          # (split outer in halves)
      [nw, xyzw, d1] = \
          slice2 (- pyz, nw, xyzw, None)

      elapsed = timer_ (elapsed)
      timer_print ("slicing time", elapsed - elapsed0)

      fma ()
      print "split_palette,\"earth.gp\" -- generate palette for pl3tree"
      split_palette ("earth.gp")
      print "gnomon -- turn on gnomon"
      gnomon (1)

      print "pl3tree with 1 slicing plane, 2 isosurfaces"
      clear3 ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      pl3tree (np, xyzp, vp, pyz)
      pl3tree (nvb, xyzvb)
      pl3tree (nwb, xyzwb)
      pl3tree (nv, xyzv)
      pl3tree (nw, xyzw)
      orient3 ()
      light3 (diffuse = .2, specular = 1)
      limits ()
      limits (square=1)
      demo5_light (1)
      paws ()
      hcp ()

      print "spin3 animated rotation, use rot3 or orient3 for one frame"
      # don't want limits to autoscale during animation
      lims = limits ( )
      spin3 ()
      limits ( ) # back to autoscaling
      demo5_light (1)
      paws ()

      light3 ()
      gnomon (0)
      limits (square = 1)
      palette ("gray.gp")

   if len (itest) == 0 or itest [0] == 4 :
      f = PR ('./bills_plot')
      n_nodes = f.NumNodes
      n_z = f.NodesOnZones
      x = f.XNodeCoords
      y = f.YNodeCoords
      z = f.ZNodeCoords
      c = f.ZNodeVelocity
      n_zones = f.NumZones
      # Put vertices in right order for Gist
      n_z = transpose (
         take (transpose (n_z), array ( [0, 4, 3, 7, 1, 5, 2, 6])))
      m3 = mesh3 (x, y, z, funcs = [c], verts = n_z ) # [0:10])
      [nv, xyzv, cv] = slice3 (m3, 1, None, None, 1, value = .9 * max (c) )
      pyz = plane3 ( array ([1, 0, 0], Float ), zeros (3, Float))
      pxz = plane3 ( array ([0, 1, 0], Float ), zeros (3, Float))

      # draw a colored plane first
      fma ()
      clear3 ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      [np, xyzp, vp] = slice3 (m3, pyz, None, None, 1)
      pl3tree (np, xyzp, vp, pyz, split = 0)
      palette ("rainbow.gp")
      orient3 ()
      demo5_light (1)
      paws ()
      

#     [nv, xyzv, d1] = \
#         slice2 (- pyz, nv, xyzv, None)
      [nw, xyzw, cw] = slice3 (m3, 1, None, None, 1, value = .9 * min (c) )
#     [nw, xyzw, d1] = \
#         slice2 (- pyz, nw, xyzw, None)
      [nvi, xyzvi, cvi] = slice3 (m3, 1, None, None, 1, value = .5 * min (c) )
      [nvi, xyzvi, cvi] = \
          slice2 (- pyz, nvi, xyzvi, cvi)
      [nvj, xyzvj, cvj] = slice3 (m3, 1, None, None, 1, value = .5 * max (c) )
      [nvj, xyzvj, cvj] = \
          slice2 (- pyz, nvj, xyzvj, cvj)

      fma ()
      print "gnomon -- turn on gnomon"
      gnomon (1)
      clear3 ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      pl3tree (nv, xyzv) # , cv)
      pl3tree (nw, xyzw) # , cw)
      pl3tree (nvi, xyzvi) # , cvi)
      pl3tree (nvj, xyzvj) # , cvi)
      orient3 ()
      light3 (ambient = 0, diffuse = .5, specular = 1, sdir = [0, 0, -1])
      limits (square=1)
      palette ("gray.gp")
      demo5_light (1)
      paws ()
      
      print "spin3 animated rotation, use rot3 or orient3 for one frame"
      # don't want limits to autoscale during animation
      spin3 ()
      limits ( ) # back to autoscaling
      demo5_light (1)
      paws ()

      light3 ()
      gnomon (0)
      palette ("gray.gp")

      draw3 ( 1 )
      paws ()
      clear3 ()
      del nv, xyzv, cv, nw, xyzw, cw, nvi, xyzvi, cvi, nvj, xyzvj, cvj
      # Make sure we don't draw till ready
      set_draw3_ (0)
      for i in range (8) :
         [nv, xyzv, cv] = slice3 (m3, 1, None, None, 1, value = .9 * min (c) +
             i * (.9 * max (c) - .9 * min (c)) / 8.)
         [nv, xyzv, d1] = \
             slice2 (pxz, nv, xyzv, None)
         pl3tree (nv, xyzv)
      orient3 ()
      light3 (ambient = 0, diffuse = .5, specular = 1, sdir = [0, 0, -1])
      limits (square=1)
      palette ("heat.gp")
      demo5_light (1)
      paws ()
      spin3 ()
      limits ( ) # back to autoscaling
      demo5_light (1)
      paws ()
      demo5_light (1)
      paws ()

   if len (itest) == 0 or itest [0] == 5 :
      # Try bert's data
      f = PR ('./berts_plot')
      nums = array ( [63, 63, 49], Int)
      dxs = array ( [2.5, 2.5, 10.], Float )
      x0s = array ( [-80., -80., 0.0], Float )
      c = f.c

      m3 = mesh3 (nums, dxs, x0s, funcs = [transpose (c)])
      [nv, xyzv, dum] = slice3 (m3, 1, None, None, value = 6.5)
      fma ()
      clear3 ()
      print "gnomon -- turn on gnomon"
      gnomon (1)
      # Make sure we don't draw till ready
      set_draw3_ (0)
      palette ("rainbow.gp")
      pl3tree (nv, xyzv)
      orient3 ()
      light3 (diffuse = .2, specular = 1)
      limits (square=1)
      demo5_light (1)
      paws ()
      spin3 ()
      demo5_light (1)
      paws ()
   if len (itest) == 0 or itest [0] == 6 :
      # Try Bill's irregular mesh
      f = PR ("ball.s0001")
      ZLss = f.ZLstruct_shapesize
      ZLsc = f.ZLstruct_shapecnt
      ZLsn = f.ZLstruct_nodelist
      x = f.sap_mesh_coord0
      y = f.sap_mesh_coord1
      z = f.sap_mesh_coord2
      c = f.W_vel_data
      # Now we need to convert this information to avs-style data
      istart = 0 # beginning index into ZLstruct_nodelist
      NodeError = "NodeError"
      ntet = 0
      nhex = 0
      npyr = 0
      nprism = 0
      nz_tet = []
      nz_hex = []
      nz_pyr = []
      nz_prism = []
      for i in range (4) :
         if ZLss [i] == 4 : # TETRAHEDRON
            nz_tet = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
                     (ZLsc [i], ZLss [i]))
            ntet = ZLsc [i]
            istart = istart + ZLss [i] * ZLsc [i]
         elif ZLss[i] == 5 : # PYRAMID
            nz_pyr = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
                     (ZLsc [i], ZLss [i]))
            npyr = ZLsc [i]
            # Now reorder the points (bill has the apex last instead of first)
            nz_pyr = transpose (
               take (transpose (nz_pyr), array ( [4, 0, 1, 2, 3])))
            istart = istart + ZLss [i] * ZLsc [i]
         elif ZLss[i] == 6 : # PRISM
            nz_prism = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
                     (ZLsc [i], ZLss [i]))
            nprism = ZLsc [i]
            # now reorder the points (bill goes around a square face 
            # instead of traversing the opposite sides in the same direction.
            nz_prism = transpose (
               take (transpose (nz_prism), array ( [0, 1, 3, 2, 4, 5])))
            istart = istart + ZLss [i] * ZLsc [i]
         elif ZLss[i] == 8 : # HEXAHEDRON
            nz_hex = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
                     (ZLsc [i], ZLss [i]))
            # now reorder the points (bill goes around a square face 
            # instead of traversing the opposite sides in the same direction.
            nz_hex = transpose (
               take (transpose (nz_hex), array ( [0, 1, 3, 2, 4, 5, 7, 6])))
            nhex = ZLsc [i]
            istart = istart + ZLss [i] * ZLsc [i]
         else :
            raise NodeError, `ZLss[i]` + "is an incorrect number of nodes."

      m3 = mesh3 (x, y, z, funcs = [c], verts = [nz_tet, nz_pyr, nz_prism,
         nz_hex])
      [nv, xyzv, cv] = slice3 (m3, 1, None, None, 1, value = .9 * max (c) )
      pyz = plane3 ( array ([1, 0, 0], Float ), zeros (3, Float))
      pxz = plane3 ( array ([0, 1, 0], Float ), zeros (3, Float))

      # draw a colored plane first
      fma ()
      clear3 ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      [np, xyzp, vp] = slice3 (m3, pyz, None, None, 1)
      pl3tree (np, xyzp, vp, pyz, split = 0)
      palette ("rainbow.gp")
      orient3 ()
      limits (square=1)
      demo5_light (1)
      paws ()

      [nw, xyzw, cw] = slice3 (m3, 1, None, None, 1, value = .9 * min (c) )
      [nvi, xyzvi, cvi] = slice3 (m3, 1, None, None, 1, value = .1 * min (c) )
      [nvi, xyzvi, cvi] = \
          slice2 (- pyz, nvi, xyzvi, cvi)
      [nvj, xyzvj, cvj] = slice3 (m3, 1, None, None, 1, value = .1 * max (c) )
      [nvj, xyzvj, cvj] = \
          slice2 (- pyz, nvj, xyzvj, cvj)
      [nvii, xyzvii, cvii] = slice3 (m3, 1, None, None, 1,
         value = 1.e-12 * min (c) )
      [nvii, xyzvii, cvii] = \
          slice2 (- pyz, nvii, xyzvii, cvii)
      [nvjj, xyzvjj, cvjj] = slice3 (m3, 1, None, None, 1,
         value = 1.e-12 * max (c) )
      [nvjj, xyzvjj, cvjj] = \
          slice2 (- pyz, nvjj, xyzvjj, cvjj)

      fma ()
      print "gnomon -- turn on gnomon"
      gnomon (1)
      clear3 ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      pl3tree (nv, xyzv) # , cv)
      pl3tree (nw, xyzw) # , cw)
      pl3tree (nvi, xyzvi) # , cvi)
      pl3tree (nvj, xyzvj) # , cvj)
      pl3tree (nvii, xyzvii) # , cvii)
      pl3tree (nvjj, xyzvjj) # , cvjj)
      orient3 ()
      light3 (ambient = 0, diffuse = .5, specular = 1, sdir = [0, 0, -1])
      limits (square=1)
      palette ("gray.gp")
      demo5_light (1)
      paws ()
      palette ("heat.gp")
      paws ()


   if len (itest) == 0 or itest [0] == 7 :
      # test plwf on the sombrero function
      # compute sombrero function
      x = arange (-20, 21, typecode = Float)
      y = arange (-20, 21, typecode = Float)
      z = zeros ( (41, 41), Float)
      r = sqrt (add.outer ( x ** 2, y **2)) + 1e-6
      z = sin (r) / r
      fma ()
      clear3 ()
      gnomon (0)
      # Make sure we don't draw till ready
      set_draw3_ (0)
      palette ("rainbow.gp")
      limits (square=1)
      orient3 ()
      light3 ()
      plwf (z, fill = z, ecolor = "black")
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      ##### Try smooth contours, log mode
      [nv, xyzv, dum] = slice3mesh (x, y, z)
      zmult = max (max (abs (x)), max (abs (y)))
      plzcont (nv, xyzv, contours = 20, scale = "normal")
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      plzcont (nv, xyzv, contours = 20, scale = "lin", edges=1)
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      plwf (z, fill = z, shade = 1, ecolor = "black")
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      plwf (z, fill = z, shade = 1, edges = 0)
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      light3(diffuse=.2,specular=1)
      print "light3,sdir=[cos(theta),.25,sin(theta)]  -- movie"
      making_movie = 1
      movie(demo5_light, lims = [xmin, xmax, ymin, ymax])
      making_movie = 0
      fma()
      demo5_light(1)
      paws ()
      plwf (z, fill = None, shade = 1, edges = 0)
      [xmin, xmax, ymin, ymax] = draw3 (1)
      palette("gray.gp")
      limits (xmin, xmax, ymin, ymax)
      paws ()


   if len (itest) == 0 or itest [0] == 8 :
      # test pl3surf on the sombrero function
      # compute sombrero function
      nc1 = 100
      nv1 = nc1 + 1
      br = - (nc1 / 2)
      tr = nc1 / 2 + 1
      x = arange (br, tr, typecode = Float) * 40. / nc1
      y = arange (br, tr, typecode = Float) * 40. / nc1
      z = zeros ( (nv1, nv1), Float)
      r = sqrt (add.outer ( x ** 2, y **2)) + 1e-6
      z = sin (r) / r
      # In order to use pl3surf, we need to construct a mesh
      # using mesh3. The way I am going to do that is to define
      # a function on the 3d mesh so that the sombrero function
      # is its 0-isosurface.
      z0 = min (ravel (z))
      z0 = z0 - .05 * abs (z0)
      maxz = max (ravel (z))
      maxz = maxz + .05 * abs (maxz)
      zmult = max (max (abs (x)), max (abs (y)))
      dz = (maxz - z0)
      nxnynz = array ( [nc1, nc1, 1], Int)
      dxdydz = array ( [1.0, 1.0, zmult*dz], Float )
      x0y0z0 = array ( [float (br), float (br), z0*zmult], Float )
      meshf = zeros ( (nv1, nv1, 2), Float )
      meshf [:, :, 0] = zmult*z - (x0y0z0 [2])
      meshf [:, :, 1] = zmult*z - (x0y0z0 [2] + dxdydz [2])

      m3 = mesh3 (nxnynz, dxdydz, x0y0z0, funcs = [meshf])
      fma ()
      # Make sure we don't draw till ready
      set_draw3_ (0)
      pldefault(edges=0)
      [nv, xyzv, col] = slice3 (m3, 1, None, None, value = 0.)
      orient3 ()
      pl3surf (nv, xyzv)
      lim = draw3 (1)
      limits (lim [0], lim [1], 1.5*lim [2], 1.5*lim [3])
      palette ("gray.gp")
      paws ()
      # Try new slicing function to get color graph
      [nv, xyzv, col] = slice3mesh (nxnynz [0:2], dxdydz [0:2], x0y0z0 [0:2],
         zmult * z, color = zmult * z)
      pl3surf (nv, xyzv, values = col)
      lim = draw3 (1)
      dif = 0.5 * (lim [3] - lim [2])
      limits (lim [0], lim [1], lim [2] - dif, lim [3] + dif)

      palette ("rainbow.gp")
      paws ()
      palette ("heat.gp")
      # Try plzcont--see if smooth mode possible
      plzcont (nv, xyzv)
      draw3 (1)
      paws ()
      plzcont (nv, xyzv, contours = 20)
      draw3 (1)
      paws ()
      plzcont (nv, xyzv, contours = 20, scale = "log")
      draw3(1)
      paws ()
      plzcont (nv, xyzv, contours = 20, scale = "normal")
      draw3(1)
      paws ()
   if len (itest) == 0 or itest [0] == 9 :
      vsf = 0.
      c = 1
      s = 1000.
      kmax = 25
      lmax = 35
      # The following computations define an interesting 3d surface.

      xr = multiply.outer (
         arange (1, kmax + 1, typecode = Float), ones (lmax, Float))
      yr = multiply.outer (
         ones (kmax, Float), arange (1, lmax + 1, typecode = Float))
      zt = 5. + xr + .2 * rand (kmax, lmax)   # ranf (xr)
      rt = 100. + yr + .2 * rand (kmax, lmax)   # ranf (yr)
      z = s * (rt + zt)
      z = z + .02 * z * rand (kmax, lmax)   # ranf (z)
      ut = rt/sqrt (rt ** 2 + zt ** 2)
      vt = zt/sqrt (rt ** 2 + zt ** 2)
      ireg =  multiply.outer ( ones (kmax, Float), ones (lmax, Float))
      ireg [0:1, 0:lmax]=0
      ireg [0:kmax, 0:1]=0
      ireg [1:15, 7:12]=2
      ireg [1:15, 12:lmax]=3
      ireg [3:7, 3:7]=0
      freg=ireg + .2 * (1. - rand (kmax, lmax))  # ranf (ireg))
      freg=array (freg, Float)
      #rt [4:6, 4:6] = -1.e8
      z [3:10, 3:12] = z [3:10, 3:12] * .9
      z [5, 5] = z [5, 5] * .9
      z [17:22, 15:18] = z [17:22, 15:18] * 1.2
      z [16, 16] = z [16, 16] * 1.1
      orient3 ()
      plwf (freg, shade = 1, edges = 0)
      [xmin, xmax, ymin, ymax] = draw3 (1)
      limits (xmin, xmax, ymin, ymax)
      paws ()
      nxny = array ( [kmax - 1, lmax - 1])
      x0y0 = array ( [0., 0.])
      dxdy = array ( [1., 1.])
      [nv, xyzv, col] = slice3mesh (nxny, dxdy, x0y0, freg)
      [nw, xyzw, col] = slice3mesh (nxny, dxdy, x0y0, freg + ut)
      pl3tree (nv, xyzv)
      pl3tree (nw, xyzw)
      draw3 (1)
      limits ( )
      paws ()

      light3 (ambient = 0, diffuse = .5, specular = 1, sdir = [0, 0, -1])
      demo5_light (1)
      paws ()

      [nv, xyzv, col] = slice3mesh (nxny, dxdy, x0y0, freg, color = freg)
      pl3surf (nv, xyzv, values = col)
      draw3 (1)
      palette ("rainbow.gp")
      paws ()
      [nv, xyzv, col] = slice3mesh (nxny, dxdy, x0y0, freg, color = z)
      pl3surf (nv, xyzv, values = col)
      draw3 (1)
      paws ()
      palette ("stern.gp")
      paws ()
      [nv, xyzv, col] = slice3mesh (nxny, dxdy, x0y0, z, color = z)
      pl3surf (nv, xyzv, values = col)
      orient3(phi=0,theta=0)
      draw3 (1)
      paws ()
      set_draw3_ (0)
      palette ("gray.gp")
      light3 ( diffuse=.1, specular=1., sdir=array([0,0,-1]))
      pl3surf (nv, xyzv)
      draw3 (1)
      paws ()

#     spin3 ()
#     paws ()

   hcp_finish ()
