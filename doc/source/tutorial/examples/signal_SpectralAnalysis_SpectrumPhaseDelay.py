import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from scipy import signal
from scipy.fft import rfft, rfftfreq

# Create input signal:
n = 50
x = np.zeros(n)
x[0] = n

# Apply FIR filter which delays signal by 3 samples:
y = signal.lfilter([0, 0, 0, 1], 1, x)

X, Y = (rfft(z_) / n for z_ in (x, y))
f = rfftfreq(n, 1)  # sampling interval T = 1 s

fig = plt.figure(tight_layout=True, figsize=(6., 4.))
gs = gridspec.GridSpec(3, 1)
ax0 = fig.add_subplot(gs[0, :])
ax1 = fig.add_subplot(gs[1:, :], sharex=ax0)

for Z_, n_, m_ in zip((X, Y), ("Input $X(f)$", "Output $Y(f)$"), ('+-', 'x-')):
    ax0.plot(f, abs(Z_), m_, alpha=.5, label=n_)
    ax1.plot(f, np.unwrap(np.angle(Z_)), m_, alpha=.5, label=n_)

ax0.set(title="Frequency Response of 3 Sample Delay Filter (no window)",
        ylabel="Magnitude", xlim=(0, f[-1]), ylim=(0, 1.1))
ax1.set(xlabel=rf"Frequency $f$ in Hertz ($\Delta f = {1/n}\,$Hz)",
        ylabel=r"Phase in rad")
ax1.set_yticks(-np.arange(0, 7)*np.pi/2,
               ['0', '-π/2', '-π', '-3/2π', '-2π', '-4/3π', '-3π'])

ax2 = ax1.twinx()
ax2.set(ylabel=r"Phase in Degrees", ylim=np.rad2deg(ax1.get_ylim()),
        yticks=np.arange(-540, 90, 90))
for ax_ in (ax0, ax1):
    ax_.legend()
    ax_.grid()
plt.show()
