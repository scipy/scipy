#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Example use of the maximum entropy package for a classification task.

    An extension of the machine translation example from the paper 'A maximum
    entropy approach to natural language processing' by Berger et al., 1996.

    Consider the translation of the English word 'in' into French.  Suppose we
    notice the following facts in a corpus of parallel texts:
        
        (1)    p(dans) + p(en) + p(à) + p(au cours de) + p(pendant) = 1
        (2)    p(dans | next English word = 'a' or 'the') = 8/10
        (3)    p(dans | c) + p(à | c)  = 1/2   for all other c
        
    This code finds the probability distribution with maximal entropy
    subject to these constraints.
"""

__author__ =  'Ed Schofield'
__version__=  '2.1'


import math
from scipy import maxentropy, sparse
import numpy

# a_grave = u'\u00e0'

samplespace = ['dans', 'en', 'à', 'au cours de', 'pendant']
contexts = ['happy', 'healthy', 'harold', 'a', 'the', 'beans']
# Occurrences of French words, and their 'next English word' contexts, in
# a hypothetical parallel corpus:
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
    return (x=='dans' or x=='à') and c not in ['a', 'the']

f = [f0, f1, f2]

numcontexts = len(contexts)
numsamplepoints = len(samplespace)

# # Dense array version:
# F = numpy.array([[f_i(x, c) for c in contexts for x in samplespace] for f_i in f])

# NEW: Sparse matrix version:
# Sparse matrices are only two dimensional in SciPy.  Store as m x size, where
# size is |W|*|X|.  Constructing this matrix by looping over all features,
# contexts, and points in the sample space, as here, is possibly very slow.
# This can and should be done more efficiently, but doing so requires some
# knowledge of the structure of the features -- so only the non-zeros can be
# looped over.
F = sparse.lil_matrix((len(f), numcontexts * numsamplepoints))
for i, f_i in enumerate(f):
    for c, context in enumerate(contexts):
        for x, samplepoint in enumerate(samplespace):
            # print "context: %s; \t sample point: %s" % (samplepoint, context)
            F[i, c * numsamplepoints + x] = f_i(samplepoint, context)

# OLD:
# Sparse matrix version
# # Sparse matrices are only two dimensional in SciPy.  Store as a list of sparse matrices,
# # one matrix for each feature f_i.
# Fs = []
# for i, f_i in enumerate(f):
#     F = sparse.lil_matrix((len(contexts), len(samplespace)))
#     for c, context in enumerate(contexts):
#         for x, samplepoint in samplespace):
#             F[c, x] = f_i(samplepoint, context)
#     Fs.append(F)
# This would work fine, and maybe have a simpler interface, but would be less
# efficient, since we'd need m=numfeatures=len(f) matrix products each
# iteration instead of just one.


# The indices_context parameter is not longer necessary if the features are
# required to be stored with contiguous indices.  This makes perfect sense with
# sparse matrices.
#
# OLD:
# # These can be in any order. The indices_context parameter to the
# # conditionalmodel constructor records this order, so indices_context[0] is an
# # array of indices all labels x in context w=0.  The simplest ordering is:
# #     (w0, x0), (w0, x1), ..., (w0, x{n-1}), (w1, x0), ...
# # in which case the indices into each row of F are:
# indices_context = numpy.arange(len(contexts) * len(samplespace)).reshape( \
#                                                 len(contexts), len(samplespace))

# The desired feature expectations, estimated from a corpus or training dataset, are
# constructed from the matrix F of feature values and the empirical pmf
# p_tilde.  This is an array with |w| * |x| elements.

# # Define a utility function
# def countall(seq):
#     """Returns a dictionary whose keys are the unique elements of seq and whose
#     values are the number of times each element occurs in seq.
#     """
#     d = {}
#     for element in seq:
#         d[element] = 1 + d.get(element, 0)
#     return d

# Compute the number N(x | c) of occurrences of each point x in each context in the training set
# # Dense vector version:  
# corpus_dict = countall(corpus)  # quicker lookups of counts
# N = numpy.array([corpus_dict.get((x, c), 0) for c in contexts for x in samplespace])
# # In general we'll want a sparse vector for this, because the length |X| * |C|
# # could be very large, but with mostly non-zeros.


# Sparse row-vector version: not yet supported by scipy.sparse:

# Utility data structures: store the indices of each context and label in a
# dict for fast lookups of their indices in their respective lists:
samplespace_index = dict((x, i) for i, x in enumerate(samplespace))
context_index = dict((c, i) for i, c in enumerate(contexts))

# Use a sparse matrix for the counts, of dimensions (1 x size), where size is
# |W| x |X|.  The element N[0, i*numcontexts+x] is the number of occurrences of
# x in context c in the training data.
N = sparse.lil_matrix((1, numcontexts * len(samplespace)))   # initialized to zero
for (x, c) in corpus:
    N[0, context_index[c] * numsamplepoints + samplespace_index[x]] += 1

# OLD:
# # Use a sparse matrix of size C x X, whose ith row vector contains all
# # points x_j in the sample space X in context c_i:
# N = sparse.lil_matrix((len(contexts), len(samplespace)))   # initialized to zero
# for (c, x) in corpus:
#     N[c, x] += 1

# This is a nicer input format, but it's more efficient internally as one long
# row vector.  Ideally sparse matrices would offer a .reshape method so this
# conversion can be done internally and transparently.


# Create a model
# Ideally, if any contexts were completely blank, we'd remove them from the
# matrix F etc first ...
model = maxentropy.conditionalmodel(F, N, numcontexts)

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
print p[c*numsamplepoints:(c+1)*numsamplepoints]

print "\nFitted distribution is:"
print "c \ x \t",
for label in samplespace:
    print label + "\t",

for c, context in enumerate(contexts):
    print "\n" + context + "\t",
    for x, label in enumerate(samplespace):
        print ("%.3f" % p[c*numsamplepoints+x]) + "\t",

print

