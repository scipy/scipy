# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from ezplot import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

print "This program tests ezplot. Type demo () to begin the test."
paws ()

# ezcshow("false")
def demo ():
   import scipy
   rand = scipy.rand
   vsf = 0.
   c = 1
   s = 1000.
   kmax = 25
   lmax = 35
   # The following computations define an interesting 3d surface.
 
   xr = multiply.outer ( arange (1, kmax + 1), ones (lmax))
   yr = multiply.outer ( ones (kmax), arange (1, lmax + 1))
   zt = 5. + xr + .2 * rand (kmax, lmax)   # ranf (xr)
   rt = 100. + yr + .2 * rand (kmax, lmax)   # ranf (yr)
   z = s * (rt + zt)
   z = z + .02 * z * rand (kmax, lmax)   # ranf (z)
   ut = rt/sqrt (rt ** 2 + zt ** 2)
   vt = zt/sqrt (rt ** 2 + zt ** 2)
   ireg =  array(multiply.outer ( ones (kmax), ones (lmax)),'i')
   ireg [0:1, 0:lmax]=0
   ireg [0:kmax, 0:1]=0
   ireg [1:15, 7:12]=2
   ireg [1:15, 12:lmax]=3
   ireg [3:7, 3:7]=0
   freg=ireg + .2 * (1. - rand (kmax, lmax))  # ranf (ireg))
   freg=freg.astype (Float)
   #rt [4:6, 4:6] = -1.e8
   z [3:10, 3:12] = z [3:10, 3:12] * .9
   z [5, 5] = z [5, 5] * .9
   z [17:22, 15:18] = z [17:22, 15:18] * 1.2
   z [16, 16] = z [16, 16] * 1.1
   cgm("close")
   win("on")
   set_mesh(rt=rt,zt=zt,ireg=ireg)
   titleb("An interesting mesh")
   plotm()
   nf()
   paws()
   titleb("Plot with ranges")
   plotm(kstyle="dashed", lstyle = "dotted", krange = (0, 20), lrange = (0,10),
      color = "yellow", thick = 2)
   nf ( )
   paws()
   titleb("L lines only")
   plotb ( )
   plotm(kstyle="none", lstyle="dashed", color="orange")
   nf ( )
   paws ()
   k2 = 3
   l2 = 7
   plotb ( )
   titleb("Region 2 singled out")
   plotm (region = 2, kstyle = "dashed", lstyle = "dotted", color = "green")
   text("Region 2", zt [k2, l2], rt [k2, l2], 32.)
   nf ()
   paws ( )
   titlet("A different contour plot")
   plotc(pvar = z, lev = 20)
   nf ()
   paws ( )
   titlet("Same plot filled.")
   plotc(pvar = z, lev = 20, color = "fill")
   nf ()
   paws ( )
   clear_mesh ( )
   nz=20
   xs = array ( range (-nz-1, nz+1) )
   ys = xs
   zs = zeros ( (2*nz + 2, 2*nz + 2), Float)
   for i in range ( len (xs) ) :
      for j in range ( len (ys)) :
         rs = sqrt ( xs [i] * xs[i] + ys [j] * ys [j] ) + 1.e-12
         zs [i, j] = sin (rs) / rs
   titleb("Sombrero, contour plot")
   plotc (pvar = zs, zt = ys, rt = xs, lev = 12)
   nf ( )
   paws ( )
   titleb("Sombrero, logarithmic contour plot")
   plotc (pvar = zs, zt = ys, rt = xs, lev = -12)
   nf ( )
   paws ( )
   titleb("Sombrero, filled contour plot")
   plotc (pvar = zs, zt = ys, rt = xs, lev = 12, color = "fill")
   nf ( )
   paws ( )
   attr (labels = "yes")
   titleb("Sombrero, contour plot using plotz.")
   plotz (zs, xs, ys, lev = 12)
   nf ( )
   paws ( )
   attr (labels = "no")
   titleb("Sombrero, filled mesh, linear")
   plotf (zs, ys, xs)
   nf ( )
   paws ( )
   titleb("Sombrero, filled mesh, logarithmic")
   plotf (zs, ys, xs, cscale = "log")
   nf ( )
   paws ( )
   titleb("Sombrero, filled mesh, normal")
   plotf (zs, ys, xs, cscale = "normal")
   nf ( )
   paws ()
   nz = 80
   xs = array ( range (-nz-1, nz+1) )/4.
   ys = xs
   zs = zeros ( (2*nz + 2, 2*nz + 2), Float)
   for i in range ( len (xs) ) :
      for j in range ( len (ys)) :
         rs = sqrt ( xs [i] * xs[i] + ys [j] * ys [j] ) + 1.e-12
         zs [i, j] = sin (rs) / rs
   titleb("Sombrero, cell array plot")
   ploti (zs)
   nf ( )
   paws ( )
   # Vector stuff
   x = zeros ( (10, 10), Float)
   y = zeros ( (10, 10), Float)
   vx = zeros ( (10, 10), Float)
   vy = zeros ( (10, 10), Float)
   ireg = zeros ( (10, 10), Int)
   for i in range (10) :
      for j in range (10) :
         vy [i, j] = sin (i)
         vx [i, j] = cos (j)
         y [i, j] = j
         x [i, j] = i
   ireg [1:5, 1:5] = 1
   ireg [1:5, 5:10] = 2
   ireg [5:10, 1:5] = 3
   ireg [5:10, 5:10] = 4
   titleb ("A vector plot")
   plotv (y, x, vy, vx, ireg, vsc = 0.5, thick = 0.1, color = "red")
   nf ( )
   paws ( )
   set_mesh (rt = y, zt = x, ut = vy, vt = vx, ireg = ireg)
   titleb ("Vector plot, regions 1 and 4.")
   plotv (vsc = 0.35, region = [1,4], color = "blue")
   nf ( )
   paws ( )
