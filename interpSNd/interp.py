import interpolateSNd as SN
import numpy as np

reload(SN)

points = np.array([[ 1.,0.,1.,0.],[0.,0.,1.,1.]])
z = np.array([1.,0.,2.,1.])
interp = SN.InterpolateSNd(points,z)

print "triangulation:\n",interp._triangulation

X=np.array([[1.4,.1,.55],[.4,.1,.3]])
print "X values:\n",X
# last component is .1 too much

output = interp(X)

print "output:\n", output