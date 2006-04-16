#!/usr/bin/env python

import scipy.special.cephes as cephes

val = 0.8
val2 = 0.3
names = dir(cephes)
for k in names[4:]:
    exec("fun = cephes."+k)
    try:
        print k + "(%f,%f) = " % (val2,val)
        print "  %f" % fun(val2,val)
    except:
        try:
            print k + "(%f) = " % val
            print "  %f" % fun(val)
        except:
            print "Error: " + k
