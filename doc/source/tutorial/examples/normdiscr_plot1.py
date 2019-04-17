import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

npoints = 20  # number of integer support points of the distribution minus 1
npointsh = npoints // 2
npointsf = float(npoints)
nbound = 4  # bounds for the truncated normal
normbound = (1 + 1/npointsf) * nbound  # actual bounds of truncated normal
grid = np.arange(-npointsh, npointsh+2, 1)  # integer grid
gridlimitsnorm = (grid-0.5) / npointsh * nbound  # bin limits for the truncnorm
gridlimits = grid - 0.5
grid = grid[:-1]
probs = np.diff(stats.truncnorm.cdf(gridlimitsnorm, -normbound, normbound))
gridint = grid
normdiscrete = stats.rv_discrete(
                        values=(gridint, np.round(probs, decimals=7)),
                        name='normdiscrete')

n_sample = 500
np.random.seed(87655678)  # fix the seed for replicability
rvs = normdiscrete.rvs(size=n_sample)
f, l = np.histogram(rvs, bins=gridlimits)
sfreq = np.vstack([gridint, f, probs*n_sample]).T
fs = sfreq[:,1] / float(n_sample)
ft = sfreq[:,2] / float(n_sample)
nd_std = np.sqrt(normdiscrete.stats(moments='v'))

ind = gridint  # the x locations for the groups
width = 0.35       # the width of the bars

plt.subplot(111)
rects1 = plt.bar(ind, ft, width, color='b')
rects2 = plt.bar(ind+width, fs, width, color='r')
normline = plt.plot(ind+width/2.0, stats.norm.pdf(ind, scale=nd_std),
                    color='b')

plt.ylabel('Frequency')
plt.title('Frequency and Probability of normdiscrete')
plt.xticks(ind+width, ind)
plt.legend((rects1[0], rects2[0]), ('true', 'sample'))

plt.show()
