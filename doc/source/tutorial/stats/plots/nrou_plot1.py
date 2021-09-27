import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sc
import math


def u_bound(x, p, center):
    if x < 0:
        return 0
    return (x - center) * x**((p-1)/2) * math.exp(-x/2)


# bounding rectangle for Gamma(p) shifted by center
def rectangle(p, center):
    h = (p+1+center)/2
    k = np.sqrt(h**2 - center*(p-1))
    u_min, u_max = u_bound(h-k, p, center), u_bound(h+k, p, center)
    v_max = math.sqrt((p-1)**(p-1) * math.exp(-(p-1)))
    return u_min, u_max, v_max


ps = np.arange(1.0, 2.5, 0.1)
reject_const, reject_const_shift = [], []
for p in ps:
    # no shift (center=0)
    u_min, u_max, v_max = rectangle(p, 0)
    reject_const.append(2*v_max*(u_max - u_min) / sc.gamma(p))
    # mode shift (center=p-1)
    u_min, u_max, v_max = rectangle(p, p-1)
    reject_const_shift.append(2*v_max*(u_max - u_min) / sc.gamma(p))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ps, reject_const, 'o', label='center=0')
ax.plot(ps, reject_const_shift, 'o', label='center=mode')
plt.xlabel('Shape parameter p of the Gamma distribution')
plt.title('NROU rejection constants for the Gamma distribution')
plt.legend()
plt.show()
