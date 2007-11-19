from scipy import *
from scipy.sandbox import montecarlo

v = rand(10)

s = montecarlo.intsampler(v)
for i in range(10):
    temp = s.sample(10**7)
