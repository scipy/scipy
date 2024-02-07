import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import rfft, rfftfreq

n, T = 100, 0.01  # number of samples and sampling interval
fcc = (20, 20.5)  # frequencies of sines
t = np.arange(n) * T
xx = (np.sin(2 * np.pi * fx_ * t) for fx_ in fcc)  # sine signals

f = rfftfreq(n, T)  # frequency bins range from 0 Hz to Nyquist freq.
XX = (rfft(x_) / n for x_ in xx)  # one-sided magnitude spectrum

fg1, ax1 = plt.subplots(1, 1, tight_layout=True, figsize=(6., 3.))
ax1.set(title=r"Magnitude Spectrum (no window) of $x(t) = \sin(2\pi f_x t)$ ",
        xlabel=rf"Frequency $f$ in Hertz (bin width $\Delta f = {f[1]}\,$Hz)",
        ylabel=r"Magnitude $|X(f)|/\tau$", xlim=(f[0], f[-1]))
for X_, fc_, m_ in zip(XX, fcc, ('x-', '.-')):
    ax1.plot(f, abs(X_), m_, label=rf"$f_x={fc_}\,$Hz")

ax1.grid(True)
ax1.legend()
plt.show()
