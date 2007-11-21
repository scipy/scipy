from scipy import *
from scipy.sandbox import montecarlo

k = [str(x) for x in range(10**6)]
v = rand(10**6)
d = dict(zip(k, v))

s = montecarlo.dictsampler(d)
for i in range(10):
    temp = s.sample(10**7)
