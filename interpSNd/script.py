# example script

import dewall as dw
import numpy as np
from numpy import array
import matplotlib.pyplot as pyplot

Pa=[array([1.1,1.1]), array([-1.,1.]), array([-1.,-1.]), \
        array([1.,-1.]), array([1.5,1.5])]

def segments(T):
    seg0 = [ [T[0][0],T[1][0]] , [T[0][1],T[1][1]] ]
    seg1 = [ [T[0][0],T[2][0]] , [T[0][1],T[2][1]] ]
    seg2 = [ [T[1][0],T[2][0]] , [T[1][1],T[2][1]] ]
    return [seg0, seg1, seg2]
    
def plot_circle(circle):
    center, radius = circle
    t = np.linspace(0,2*np.pi,500)
    x_offset = radius*np.sin(t)
    y_offset = radius*np.cos(t)
    x = center[0] + x_offset
    y = center[1] + y_offset
    pyplot.plot(x,y)
    
Pb=[]#[array([.25, -.25]), array([0,.75])]

P = Pa+Pb

P = [ np.array([np.random.gamma(1), np.random.gamma(1)]) \
        for j in range(10) ]

triangul = dw.dewall(P)

# plotting the known data points
for p in P: pyplot.scatter([p[0]],[p[1]])

# plotting the triangulation
for tri in triangul:
    for seg in segments(tri):
        pyplot.plot(seg[0],seg[1])
        
# plotting the circumcircles
for tri in triangul:
    pass#plot_circle(dw.circumcircle(tri))
        
pyplot.show()

print "triangulation:\n",triangul