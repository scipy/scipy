#!/usr/bin/env python
#  $Id$
#  ---------------------------------------------------------------------
#
#  NAME:     gistdemomovie.py
#
#  PURPOSE:  Mesh plotting demo
#  Adapted from demo2.i for yorick by Michiel de Hoon
#
#  CHANGES: 
#  05/05/03 mdh Original conversion. 
#
#  ---------------------------------------------------------------------
#    Copyright (c) 1995.  The Regents of the University of California.
#                   All rights reserved.

from gist import *

def run(which=None, time_limit=60):
   """Exhibit quadrilateral mesh plots in 3 movies of a drumhead.
     The drumhead is initially stationary, but has a bump near one
     edge.  Yorick is solving a 2D wave equation to compute the
     evolution of this bump.

     The first movie is a filled mesh plot with color "proportional"
     to the height of the surface of the drum.  A few well chosen
     contour levels (here 3) add a lot to a filled mesh plot.

     The second movie is a "3D" perspective plot of the height of the
     drumhead.  In this movie, the mesh lines are drawn, which is
     slightly confusing since the cells are not all the same shape.

     The second movie is a "3D" shaded plot of the height of the
     drumhead.  Yorick computes surface shading based on the angle
     of each cell from a light source.

     As you watch this, you might reflect on the two dimensionality
     of your retina.  What Yorick lacks by way of 3D graphics is
     really just fancy hidden surface algorithms; the simple
     painter's algorithm used here and in plwf.py is easy to
     implement.

     There are two optional arguments to demo2: the first is the
     number of the movie (1, 2, or 3) you want to watch; the second
     is a time limit on the duration of each movie in seconds (default
     is 60 seconds each)."""
   import movie
   global f, fdot, dt, x, y, level

   # generate a 30-by-30 cell mesh on the [-1,1] square
   x= span(-1, 1, 31, 31)
   y= transpose(x)
   # map the square mesh into a mesh on the unit circle
   # this mesh has more nearly equal area cells than a polar
   # coordinate circle
   scale= maximum(abs(y),abs(x))/(hypot(y,x)+1.e-30)
   x= x*scale
   y= y*scale

   f= exp(-8.*hypot(y+.67,x+.25)**2)*(1.-hypot(y,x)**2)
   f0 = array(f) # get an independent copy
   fdot= 0.0*f[1:-1,1:-1]

   lf= laplacian(f, y,x)
   xdz= x[1:,1:]+x[:-1,1:]-x[1:,:-1]-x[:-1,:-1]
   xzd= x[1:,1:]-x[:-1,1:]+x[1:,:-1]-x[:-1,:-1]
   ydz= y[1:,1:]+y[:-1,1:]-y[1:,:-1]-y[:-1,:-1]
   yzd= y[1:,1:]-y[:-1,1:]+y[1:,:-1]-y[:-1,:-1]
   dt= 0.1875*sqrt(min(min(abs(xdz*yzd - xzd*ydz))))

   window(0, wait=1, style="nobox.gs")
   palette("heat.gp")
   limits(-1, 1, -1, 1)

   # roll the filled mesh movie
   if which==None or which==1:
     fc= (f[1:,1:]+f[:-1,1:]+f[1:,:-1]+f[:-1,:-1]) / 4.
     cmin= cmax= max([max(abs(row)) for row in fc])
     cmin= -cmin
     level= cmax/4.
     display_plf(0)
     fixedlimits = limits()
     movie.movie(display_plf, time_limit, lims=fixedlimits, timing=1)
     # Note; movie_timing is a global variable in movie.py
     print movie.movie_timing[3], "frames of filled mesh drumhead completed in",
     print movie.movie_timing[2], "sec" 
     print "Rate for filled mesh is",
     print movie.movie_timing[3]/(movie.movie_timing[0]-movie.movie_timing[4]+1.0e-4), 
     print "frames/(CPU sec),",
     print movie.movie_timing[3]/(movie.movie_timing[2]-movie.movie_timing[4]+1.0e-4),
     print "frames(wall sec)"

   # roll the perspective movie */
   if which==None or which==2:
     f[:,:]= f0
     limits(-1,1,-1,1)
     display_plm(0)
     fixedlimits = limits()
     movie.movie(display_plm, time_limit, lims=fixedlimits, timing=1)
     print movie.movie_timing[3], "frames of wireframe surface drumhead completed in",
     print movie.movie_timing[2], "sec"
     print "Rate for filled mesh is",
     print movie.movie_timing[3]/(movie.movie_timing[0]-movie.movie_timing[4]+1.0e-4), 
     print "frames/(CPU sec),",
     print movie.movie_timing[3]/(movie.movie_timing[2]-movie.movie_timing[4]+1.0e-4),
     print "frames(wall sec)"

   # roll the shaded movie
   if which==None or which==3:
     f[:,:]= f0
     limits(-1,1,-1,1)
     display_pl3(0)
     fixedlimits = limits()
     movie.movie(display_pl3, time_limit, lims=fixedlimits, timing=1)
     print movie.movie_timing[3], "frames of filled surface drumhead completed in",
     print movie.movie_timing[2], "sec"
     print "Rate for filled mesh is",
     print movie.movie_timing[3]/(movie.movie_timing[0]-movie.movie_timing[4]+1.0e-4), 
     print "frames/(CPU sec),",
     print movie.movie_timing[3]/(movie.movie_timing[2]-movie.movie_timing[4]+1.0e-4),
     print "frames(wall sec)"

     fma()
     limits()

def display_plf(i):
   # display first
   global fdot,f,level
   fc= (f[1:,1:]+f[:-1,1:]+f[1:,:-1]+f[:-1,:-1]) / 4.
   cmin= cmax= max([max(abs(row)) for row in fc])
   cmin= -cmin
   plf(fc, -y, -x, cmin=cmin, cmax=cmax)
   # the 0 contour level is too noisy without some smoothing...
   ftemp = (f[1:,1:]+f[1:,:-1]+f[:-1,1:]+f[:-1,:-1])/4.
   fs = zeros(shape(f),'d')
   fs[1:,1:] = ftemp
   fs[:-1,1:] = fs[:-1,1:] + ftemp
   fs[1:,:-1] = fs[1:,:-1] + ftemp
   fs[:-1,:-1] = fs[:-1,:-1] + ftemp
   fs[:,1:-1] = fs[:,1:-1]/2.
   fs[1:-1,:] = fs[1:-1,:]/2.

   plc(fs, levs=[0.], marks=0, color="green", type="solid")
   plc(f,  levs=[level], marks=0, color="black", type="dash")
   plc(f,  levs=[-level], marks=0, color="green", type="dash")

   # then take a step forward in time
   lf= laplacian(f, y,x)
   fdot = fdot + lf*dt
   f[1:-1,1:-1] = f[1:-1,1:-1] + fdot*dt

   return i<200

def display_plm(i):
  global fdot
  # display first
  pl3d(0, f, y, x)

  # then take a step forward in time
  lf= laplacian(f, y,x)
  fdot = fdot + lf*dt
  f[1:-1,1:-1] = f[1:-1,1:-1] + fdot*dt

  return i<200

def display_pl3(i):
  global fdot
  # display first
  pl3d(1, f, y, x)

  # then take a step forward in time
  lf= laplacian(f, y,x)
  fdot = fdot + lf*dt
  f[1:-1,1:-1] = f[1:-1,1:-1] + fdot*dt

  return i<200

def laplacian(f, y,x):
   # There are many ways to form the Laplacian as a finite difference.
   # This one is nice in Yorick.
   # Start with the two median vectors across each zone.
   fdz= (f[1:,1:]-f[1:,:-1]+f[:-1,1:]-f[:-1,:-1])/2.
   fzd= (f[1:,1:]+f[1:,:-1]-f[:-1,1:]-f[:-1,:-1])/2.
   xdz= (x[1:,1:]-x[1:,:-1]+x[:-1,1:]-x[:-1,:-1])/2.
   xzd= (x[1:,1:]+x[1:,:-1]-x[:-1,1:]-x[:-1,:-1])/2.
   ydz= (y[1:,1:]-y[1:,:-1]+y[:-1,1:]-y[:-1,:-1])/2.
   yzd= (y[1:,1:]+y[1:,:-1]-y[:-1,1:]-y[:-1,:-1])/2.

   # Estimate the gradient at the center of each cell.
   area= xdz*yzd - xzd*ydz
   gradfx= (fdz*yzd - fzd*ydz)/area
   gradfy= (xdz*fzd - xzd*fdz)/area

   # Now consider the mesh formed by the center points of the original.
   x= (x[1:,1:]+x[:-1,1:]+x[1:,:-1]+x[:-1,:-1]) / 4.
   y= (y[1:,1:]+y[:-1,1:]+y[1:,:-1]+y[:-1,:-1]) / 4.
   xdz= x[:,1:] - x[:,:-1]
   xzd= x[1:,:] - x[:-1,:]
   ydz= y[:,1:] - y[:,:-1]
   yzd= y[1:,:] - y[:-1,:]
   area= ((xdz[1:,:]+xdz[:-1,:])*(yzd[:,1:]+yzd[:,:-1]) -(xzd[:,1:]+xzd[:,:-1])*(ydz[1:,:]+ydz[:-1,:]))/4.
   term1 = xdz*(gradfy[:,1:]+gradfy[:,:-1])-ydz*(gradfx[:,1:]+gradfx[:,:-1])
   term2 = yzd*(gradfx[1:,:]+gradfx[:-1,:])-xzd*(gradfy[1:,:]+gradfy[:-1,:])

   return (term1[1:,:]-term1[:-1,:]+term2[:,1:]-term2[:,:-1]) / (2.*area)


def pl3d(shading, z, y, x):
  # rotate so that (zp,yp) are screen (y,x)
  # These orientations are cunningly chosen so that the painter's
  #  algorithm correctly draws hidden surfaces first -- see help, plf
  #  for a description of the order cells are drawn by plf.
  theta= 30. * pi/180.  # angle of viewer above drumhead
  phi= 120. * pi/180.

  ct= cos(phi)
  st= sin(phi)
  yp= y*ct - x*st
  xp= x*ct + y*st

  ct= cos(theta)
  st= sin(theta)
  zp= z*ct - xp*st
  xp= xp*ct + z*st

  if not shading:
    color= []
    edges= 1;
  else:
    # compute the two median vectors for each cell
    m0x= (xp[1:,1:]+xp[:-1,1:]-xp[1:,:-1]-xp[:-1,:-1])/2
    m0y= (yp[1:,1:]+yp[:-1,1:]-yp[1:,:-1]-yp[:-1,:-1])/2
    m0z= (zp[1:,1:]+zp[:-1,1:]-zp[1:,:-1]-zp[:-1,:-1])/2
    m1x= (xp[1:,1:]-xp[:-1,1:]+xp[1:,:-1]-xp[:-1,:-1])/2
    m1y= (yp[1:,1:]-yp[:-1,1:]+yp[1:,:-1]-yp[:-1,:-1])/2
    m1z= (zp[1:,1:]-zp[:-1,1:]+zp[1:,:-1]-zp[:-1,:-1])/2
    # define the normal vector to be their cross product
    nx= m0y*m1z - m0z*m1y
    ny= m0z*m1x - m0x*m1z
    nz= m0x*m1y - m0y*m1x
    n= sqrt(nx*nx+ny*ny+nz*nz)
    nx= nx / n
    ny= ny / n
    nz= nz / n 
    color= bytscl(nx, cmin=0.0, cmax=1.0)
    edges= 0

  plf(color, zp, yp, edges=edges)
