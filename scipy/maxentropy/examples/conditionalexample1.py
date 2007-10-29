# Test for conditional models
# Ed Schofield, 2006

from numpy import *
from scipy.maxentropy import *

# Two contexts W, four labels x
# E_p f_0(W, X) = 0.4
# where f_0(w, x) = indicator func "is the label x=0 and is the context w=0?"
# So we want the distribution:
# x \ w         0       1
# 0             0.4     0.25
# 1             0.2     0.25
# 2             0.2     0.25
# 3             0.2     0.25

# We can achieve this by creating a feature matrix with one row per constraint,
# as follows:
F = array([[1, 0, 0, 0, 0, 0, 0, 0]])
# Each column represents one (w, x) pair.  The number of columns is the product
# |w| * |x|, here 8.  The order is (w0,x0), (w0,x1), (w0, x2), ..., (w1, x0),
# etc.
numcontexts = 2
numlabels = 4

# OLD:
# These can be in any order. The indices_context parameter to the
# conditionalmodel constructor records this order, so indices_context[0] is an
# array of indices all labels x in context w=0.  The simplest ordering is:
#     (w0, x0), (w0, x1), ..., (w0, x{n-1}), (w1, x0), ...
# in which case the parameter is:
# indices_context = array([[0, 1, 2, 3], [4, 5, 6, 7]])

# The counts of each (w, x) pair, estimated from a corpus or training dataset, is
# stored as an array with |w| * |x| elements in same order as before.
counts = array([4, 3, 2, 1, 4, 3, 2, 1])
# Note that, since the only active feature was for the first element (w0, x0),
# only the first value is relevant.  The others are subject to no constraints,
# and will be chosen to maximize entropy.

model = conditionalmodel(F, counts, numcontexts)
model.verbose = True
model.fit()
# Do it again, since the CG algorithm gets stuck sometimes.  Why is this??
model.fit()
# Note: to use the bound-constrained limited memory BFGS algorithm instead, we
# would use:
# model.fit(algorithm='LBFGSB')

# Display the fitted model
pmf = model.pmf()
# The elements of this are flatted like the rows of F and p_tilde.  We display
# them nicely:
print "x \ w \t 0 \t 1",
for x in range(4):
    print '\n' + str(x),
    for w in range(2):
        print ' \t %.3f' % pmf[w*numlabels + x],
        # print ' \t %.3f' % pmf[indices_context[w]][x],
print
