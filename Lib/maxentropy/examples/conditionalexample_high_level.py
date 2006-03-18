#!/usr/bin/env python

""" THIS EXAMPLE DOESN'T YET WORK!

    Example use of the maximum entropy package for a classification task.

    An extension of the machine translation example from the paper 'A maximum
    entropy approach to natural language processing' by Berger et al., 1996.

    Consider the translation of the English word 'in' into French.  Suppose we
    notice the following facts in a corpus of parallel texts:
        
        (1)    p(dans) + p(en) + p(a) + p(au cours de) + p(pendant) = 1
        (2)    p(dans | next English word = 'a' or 'the') = 8/10
        (3)    p(dans | c) + p(a | c)  = 1/2   for all other c
        
    This code finds the probability distribution with maximal entropy
    subject to these constraints.
"""

__author__ =  'Ed Schofield'
__version__=  '2.1'


import math
from scipy import maxentropy
import numpy

a_grave = u'\u00e0'

samplespace = ['dans', 'en', a_grave, 'au cours de', 'pendant']
contexts = ['happy', 'healthy', 'harold', 'a', 'the', 'beans']
# Occurrences of French words, and their 'next English word' contexts, in a parallel corpus:
corpus = [('dans', 'a'), ('dans', 'a'), ('dans', 'a'), ('dans', 'the'), ('pendant', 'a')] + \
         [('dans', 'happy'), ('au cours de', 'healthy')]

def f0(x, c):
    return x in samplespace

def f1(x, c):
    if x == 'dans' and c in ['a', 'the']:
        return True
    else:
        return False

def f2(x, c):
    return (x=='dans' or x==a_grave) and c not in ['a', 'the']

f = [f0, f1, f2]

F = numpy.array([[f_i(x, c) for c in contexts for x in samplespace] for f_i in f])
# These can be in any order. The indices_context parameter to the
# conditionalmodel constructor records this order, so indices_context[0] is an
# array of indices all labels x in context w=0.  The simplest ordering is:
#     (w0, x0), (w0, x1), ..., (w0, x{n-1}), (w1, x0), ...
# in which case the indices into each row of F are:
indices_context = numpy.arange(len(contexts) * len(samplespace)).reshape( \
                                                len(contexts), len(samplespace))

# The desired feature expectations, estimated from a corpus or training dataset, are
# constructed from the matrix F of feature values and the empirical pmf
# p_tilde.  This is an array with |w| * |x| elements.

# Define a utility function
def countall(seq):
    """Returns a dictionary whose keys are the unique elements of seq and whose
    values are the number of times each element occurs in seq.
    """
    d = {}
    for element in seq:
        d[element] = 1 + d.get(element, 0)
    return d

corpus_dict = countall(corpus)  # quicker lookups of counts

# Compute the number N(x | c) of occurrences of each point x in each context in the training set
N = [corpus_dict.get((x, c), 0) for c in contexts for x in samplespace]
# In general we'll want a sparse vector for this, because the length |X| * |C|
# could be very large, but with mostly non-zeros.

# Create a model
# Ideally, if any contexts were completely blank, we'd remove them from the
# matrix F etc first ...
model = maxentropy.conditionalmodel(F, indices_context, N)

# Now set the desired feature expectations

model.verbose = True

# Fit the model
model.fit()
model.fit()

# Output the distribution
print "\nFitted model parameters are:\n" + str(model.theta)

p = model.probdist()

print "\npmf table p(x | c), where c is the context 'the':"
c = contexts.index('the')
print p[indices_context[c]]

print "\nFitted distribution is:"
print "c \ x \t",
for label in samplespace:
    print label + "\t",
for c, context in enumerate(contexts):
    print "\n" + context + "\t",
    for x, label in enumerate(samplespace):
        print ("%.3f" % p[indices_context[c][x]]) + "\t",

print

