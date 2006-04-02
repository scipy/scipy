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

samplespace = ['dans', 'en', 'à', 'au cours de', 'pendant']
contexts = ['happy', 'healthy', 'harold', 'a', 'the', 'beans']
# Occurrences of French words, and their 'next English word' contexts, in
# a hypothetical parallel corpus:
corpus = [('dans', 'a'), ('dans', 'a'), ('dans', 'a'), ('dans', 'the'), \
          ('pendant', 'a'), ('dans', 'happy'), ('au cours de', 'healthy')]

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

# Utility data structures: store the indices of each context and label in a
# dict for fast lookups of their indices in their respective lists:
samplespace_index = dict((x, i) for i, x in enumerate(samplespace))
context_index = dict((c, i) for i, c in enumerate(contexts))

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

# Note: just looping over the corpus, instead of the entire sample space, yields
# a different fitted distribution.  For example, using
#
# F = sparse.lil_matrix((len(f), numcontexts * numsamplepoints))
# for i, f_i in enumerate(f):
#     for (x, c) in corpus:
#         F[i, context_index[c] * numsamplepoints + samplespace_index[x]] = \
#                 f_i(x, c)
#
# yields 
#
# c \ x   dans    en      à       au cours de     pendant
# happy   1.000   0.000   0.000   0.000           0.000
# healthy 0.000   0.000   0.000   1.000           0.000
# harold  0.200   0.200   0.200   0.200           0.200
# a       0.750   0.000   0.000   0.000           0.250
# the     1.000   0.000   0.000   0.000           0.000
# beans   0.200   0.200   0.200   0.200           0.200
#
# One difference is that the distribution for the context 'beans' is uniform. We
# wanted to impose a constraint upon it, but this context never appeared in the
# corpus. So how could this have worked? To be able to impose constraints like
# this we need to compute the features of all sample points eventually.  But
# perhaps some of this could be postponed until later (e.g. testing time), rather than
# required up-front?


# Store the counts of each (context, sample point) pair in the corpus, in a
# sparse matrix of dimensions (1 x size), where size is |W| x |X|.  The element
# N[0, i*numcontexts+x] is the number of occurrences of x in context c in the
# training data.
# (The maxentropy module infers the empirical pmf etc. from the counts N)

N = sparse.lil_matrix((1, numcontexts * len(samplespace)))   # initialized to zero
for (x, c) in corpus:
    N[0, context_index[c] * numsamplepoints + samplespace_index[x]] += 1

# OLD:
# # Use a sparse matrix of size C x X, whose ith row vector contains all
# # points x_j in the sample space X in context c_i:
# N = sparse.lil_matrix((len(contexts), len(samplespace)))   # initialized to zero
# for (c, x) in corpus:
#     N[c, x] += 1

# This would be a nicer input format, but computations are more efficient
# internally with one long row vector.  Ideally sparse matrices would
# offer a .reshape method so this conversion could be done internally and
# transparently.  Then the numcontexts argument to the conditionalmodel
# constructor could also be inferred from the matrix dimensions.

# Create a model
model = maxentropy.conditionalmodel(F, N, numcontexts)

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

