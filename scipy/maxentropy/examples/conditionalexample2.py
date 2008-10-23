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

from scipy import maxentropy, sparse

samplespace = ['dans', 'en', 'à', 'au cours de', 'pendant']
# Occurrences of French words, and their 'next English word' contexts, in
# a hypothetical parallel corpus:
corpus = [('dans', 'a'), ('dans', 'a'), ('dans', 'a'), ('dans', 'the'), \
          ('pendant', 'a'), ('dans', 'happy'), ('au cours de', 'healthy')]
contexts = list(set([c for (x, c) in corpus]))

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
# dict for fast lookups of their indices into their respective lists:
samplespace_index = dict([(x, i) for i, x in enumerate(samplespace)])
context_index = dict([(c, i) for i, c in enumerate(contexts)])

# # Dense array version:
# F = numpy.array([[f_i(x, c) for c in contexts for x in samplespace] for f_i in f])

# NEW: Sparse matrix version:
# Sparse matrices are only two dimensional in SciPy.  Store as m x size, where
# size is |W|*|X|.
F = sparse.lil_matrix((len(f), numcontexts * numsamplepoints))
for i, f_i in enumerate(f):
    for c, context in enumerate(contexts):
        for x, samplepoint in enumerate(samplespace):
            F[i, c * numsamplepoints + x] = f_i(samplepoint, context)


# Store the counts of each (context, sample point) pair in the corpus, in a
# sparse matrix of dimensions (1 x size), where size is |W| x |X|.  The element
# N[0, i*numcontexts+x] is the number of occurrences of x in context c in the
# training data.
# (The maxentropy module infers the empirical pmf etc. from the counts N)

N = sparse.lil_matrix((1, numcontexts * len(samplespace)))   # initialized to zero
for (x, c) in corpus:
    N[0, context_index[c] * numsamplepoints + samplespace_index[x]] += 1

# Ideally, this could be stored as a sparse matrix of size C x X, whose ith row
# vector contains all points x_j in the sample space X in context c_i:
# N = sparse.lil_matrix((len(contexts), len(samplespace)))   # initialized to zero
# for (c, x) in corpus:
#     N[c, x] += 1

# This would be a nicer input format, but computations are more efficient
# internally with one long row vector.  What we really need is for sparse
# matrices to offer a .reshape method so this conversion could be done
# internally and transparently.  Then the numcontexts argument to the
# conditionalmodel constructor could also be inferred from the matrix
# dimensions.

# Create a model
model = maxentropy.conditionalmodel(F, N, numcontexts)

model.verbose = True

# Fit the model
model.fit()

# Output the distribution
print "\nFitted model parameters are:\n" + str(model.params)

p = model.probdist()

print "\npmf table p(x | c), where c is the context 'the':"
c = contexts.index('the')
print p[c*numsamplepoints:(c+1)*numsamplepoints]

print "\nFitted distribution is:"
print "%12s" % ("c \ x"),
for label in samplespace:
    print "%12s" % label,

for c, context in enumerate(contexts):
    print "\n%12s" % context,
    for x, label in enumerate(samplespace):
        print ("%12.3f" % p[c*numsamplepoints+x]),

print
