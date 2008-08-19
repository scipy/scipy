# example sript
import _triang as _t
import numpy as np

print "dir(_t): ", dir(_t)
a = _t.sayHello(5 ,'foo')
print "it worked?  I won't see this"
print "a equals: ", a

arr = np.array([1.0, 2.3, 4.8])
print "other answer: ", _t.triangulate(arr)
print "now at end"