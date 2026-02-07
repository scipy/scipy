import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import rfft, rfftfreq
from scipy.signal import detrend

rng = np.random.default_rng()  # seeding with default seed for reproducibility

N, tau = 1000, 0.1 # number of samples and signal duration in seconds
T = tau / N  # sampling interval
t = np.arange(N) * T  # sample times

# Create signal:
x_sines = sum(np.sin(2 * np.pi * f_ * t) for f_ in [200, 300, 400])
x_drift = 1e3 * (1 - np.cos(2*t))  # almost linear drift
x_noise = rng.normal(scale=np.sqrt(5), size=t.shape)  # Gaussian noise
x = x_drift + x_sines + x_noise  # the signal
x_d = detrend(x, type='linear')  # signal with linear trend removed

f = rfftfreq(N, T)  # frequency values
aX, aX_d = abs(rfft(x) / N), abs(rfft(x_d) / N)  # Magnitude spectrum

fig = plt.figure(figsize=(6., 5.))
ax0, ax1, ax2 = [fig.add_axes((0.11, bottom, .85, .2)) for bottom in [.7, .33, .1]]
for c_, (x_, n_) in enumerate(zip((x, x_d), ("$x(t)$", "$x_d(t)$"))):
    ax0.plot(t, x_, f'C{c_}-', alpha=0.5, label=n_)
ax0.set(title="Signal $x(t)$ and detrended signal $x_d(t)$", ylabel="Amplitude",
        xlabel=rf"Time $t$ in seconds ($T={1e3*T:g}\,$ms, {N} samples)", xlim=(0, tau))
for c_, (ax_, aX_, n_) in enumerate(zip((ax1, ax2), (aX, aX_d), ("X", "X_d"))):
    ax_.plot(f, aX_, f'C{c_}-', alpha=0.5, label=f"$|{n_}(f)|$")
    ax_.set(xlim=(0, 500), ylim=(0, 1.3), ylabel="Magnitude")
ax1.set(title="Magnitude Spectrum", xticklabels=[])
ax2.set_xlabel(rf"Frequency $f$ in hertz ($Δf={f[1]:g}\,$Hz, $f_S={1e-3/T:g}\,$kHz)")
for ax_ in (ax0, ax1, ax2):
    ax_.grid(True)
    ax_.legend()

plt.show()
