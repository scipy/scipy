import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

b = [1, 1]
w, h = signal.freqz(b)

fig, ax = plt.subplots()
ax.plot(w, 20 * np.log10(abs(h)))
ax.set(title='Frequency response of filter', xlabel='Frequency [rad/sample]', ylabel='Magnitude [dB]')
ax.grid(True)

plt.show()

