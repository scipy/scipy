"""Generates 10**8 random variates (strings) drawn from a discrete distribution
on the sample space {'a', 'b', 'c'}.

Run this script with timeit to get an idea of how the speed compares with
sampling strings and/or integers over a much larger sample space.

The point of the montecarlo module's compact 5-table sampler is for the time
for simulating variates to be independent (or nearly so) of the size of the
sample space.
"""

from scipy import *
from scipy.sandbox import montecarlo

d = {'a':0.2,'b':0.3,'c':0.5}

s = montecarlo.dictsampler(d)
for i in range(10):
    temp = s.sample(10**7)
