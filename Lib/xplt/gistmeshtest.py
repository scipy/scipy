# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from region import *
from quadmesh import *
from graph2d import *
from MLab import rand
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

vsf = 0.
c = 1
s = 1000.
kmax = 25
lmax = 35
xr = multiply.outer ( arange (1, kmax + 1, typecode = Float), ones (lmax))
yr = multiply.outer ( ones (kmax), arange (1, lmax + 1, typecode = Float))
zt = 5. + xr + .2 * rand (kmax, lmax)   # ranf (xr)
rt = 100. + yr + .2 * rand (kmax, lmax)   # ranf (yr)
z = s * (rt + zt)
z = z + .02 * z * rand (kmax, lmax)   # ranf (z)
ut = rt/sqrt (rt ** 2 + zt ** 2)
vt = zt/sqrt (rt ** 2 + zt ** 2)
ireg =  multiply.outer ( ones (kmax), ones (lmax))
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
sh = shape (z)
sh1 = sh[0]
sh2 = sh[1]
x = multiply.outer (arange (25, typecode = Float), ones (35, Float))
y = multiply.outer (ones (25, Float), arange (35, typecode = Float))

def demo () :
   qm = QuadMesh (x = zt, y = rt, ireg = ireg, width = 1., color = "blue")
   gr = Graph2d (qm)
   gr.plot ()
   paws ( )
   # qm.set (width = 3.)
   # gr.plot ()
   # paws ( )
   qm.set (color = "red", width = 1.)
   gr.change_plot(xyequal = 1)
   paws ()
   qm.set (color = "yellow")
   gr.change_plot(xyequal = 0)
   paws ()
   # Now we need two meshes, one with boundaries of regions showing,
   # the other emphasizing Region 2.
   r1 = Region (number = 1, width = 1., color = "fg", boundary = 1)
   r2 = Region (number = 2, width = 1., color = "green", type = "dash")
   r3 = Region (number = 3, width = 1., color = "fg", boundary = 1)
   
   qm.set (regions = [r1, r2, r3])
   gr.change(text = "Region 2", text_pos = [0.25,0.54], text_size = 18, text_color = "red")
   gr.plot ()
   paws ()
   qm.set (vx = vt, vy = ut, color = "red", scale = 1.)
   r1.set (color = "red")
   r2.set (color = "fg", type = "solid", boundary = 1)
   r3.set (vectors = 0)
   gr.change(text = "")
   gr.plot ()
   paws ()
   qm.set(z=z,x=xr, y = yr, vx = xr + yr/5., vy = yr + xr/10., scale = .05)
   r1.set (color = "orange",width = 3.)
   r2.set (color = "red",width = 3., vectors = 0, contours = 1, type = "solid",levels=20)
   r3.set (color = "cyan",width = 3., vectors = 1)
   gr.plot ()
   paws ()
   
   qm2 = QuadMesh (z = z, y = y, x = x, color = "yellow", width = 3.,
                   levels = 20, marks = 1)
   gr.delete (1)
   gr.add (qm2)
   gr.plot ()
   paws()
   qm2 = QuadMesh (z = z - z[kmax / 2, lmax / 2] , y = yr, x = xr, type = "dash",
                   color = "fg", levels = 20, width = 2., marks = 1)
   gr.delete (1)
   gr.add (qm2)
   gr.plot ()
   paws()
   qm2.set(filled=1)
   gr.plot ()
   paws()
   qm2 = QuadMesh (z = z, y = yr, x = xr , color = "purple" , marks = 1,
                   marker = "O", width = 3.)
   gr.delete (1)
   gr.add (qm2)
   gr.plot ()
   paws()
   qm2 = QuadMesh (z = z, y = yr, x = xr , color = "cyan" , marks = 1,
       width = 3., levels = [0.,max (ravel(z)) / 4., 3.*max (ravel(z)) / 4.,
       7. *max (ravel(z)) / 8.])
   gr.delete (1)
   gr.add (qm2)
   gr.plot ()
   paws()

print "Type demo() to proceed."
