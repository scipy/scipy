#!/usr/bin/env python

""" Example use of the maximum entropy module:

    Machine translation example -- English to French -- from the paper 'A
    maximum entropy approach to natural language processing' by Berger et
    al., 1996.

    Consider the translation of the English word 'in' into French.  We
    notice in a corpus of parallel texts the following facts:

        (1)    p(dans) + p(en) + p(a) + p(au cours de) + p(pendant) = 1
        (2)    p(dans) + p(en) = 3/10
        (3)    p(dans) + p(a)  = 1/2

    This code finds the probability distribution with maximal entropy
    subject to these constraints.
"""

__author__ =  'Ed Schofield'
__version__=  '2.1'

from scipy import maxentropy

a_grave = u'\u00e0'

samplespace = ['dans', 'en', a_grave, 'au cours de', 'pendant']

def f0(x):
    return x in samplespace

def f1(x):
    return x=='dans' or x=='en'

def f2(x):
    return x=='dans' or x==a_grave

f = [f0, f1, f2]

model = maxentropy.model(f, samplespace)

# Now set the desired feature expectations
K = [1.0, 0.3, 0.5]

model.verbose = True

# Fit the model
model.fit(K)

# Output the distribution
print "\nFitted model parameters are:\n" + str(model.params)
print "\nFitted distribution is:"
p = model.probdist()
for j in range(len(model.samplespace)):
    x = model.samplespace[j]
    print ("\tx = %-15s" %(x + ":",) + " p(x) = "+str(p[j])).encode('utf-8')


# Now show how well the constraints are satisfied:
print
print "Desired constraints:"
print "\tp['dans'] + p['en'] = 0.3"
print ("\tp['dans'] + p['" + a_grave + "']  = 0.5").encode('utf-8')
print
print "Actual expectations under the fitted model:"
print "\tp['dans'] + p['en'] =", p[0] + p[1]
print ("\tp['dans'] + p['" + a_grave + "']  = " + str(p[0]+p[2])).encode('utf-8')
# (Or substitute "x.encode('latin-1')" if you have a primitive terminal.)
