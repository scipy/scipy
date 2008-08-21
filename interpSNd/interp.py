import interpolateSNd as SN
import numpy as np
import dewall as dw

reload(SN)

if False:
    points = np.array([[ 1.,0.,1.,0.],[0.,0.,1.,1.]])
    z = np.array([1.,0.,2.,1.])
    interp = SN.InterpolateSNd(points,z)

    print "triangulation:\n",interp._triangulation

    X=np.array([[1.4,.1,.55],[.4,.1,.3]])
    print "X values:\n",X

    output = interp(X)

    print "output:\n", output
    
    
if True:
    points = np.array([  [0., 0, 0, 1., 1., 1., 0., 1., .2],
                                [0., 1., 0, 0, 1., 0., 1., 1., .2],
                                [0., 0, 1., 0, 0., 1., 1., 1., .2] ])
    z = np.sum(points,axis=0).reshape(points.shape[1])
    interp = SN.InterpolateSNd(points,z)
    
    print "*"*100+'\nMADE IT'
    
    X = np.array([ [.1,.2,.1,.1],
                        [.1,.1,.2,.1],
                        [.1,.1,.1,.2] ])
                        
    output = interp(X)
    
    print "output:\n",output
    
if False:
    P = [np.random.random_sample(3) for i in range(7)]
    print "P:",P
    tri = dw.dewall(P)