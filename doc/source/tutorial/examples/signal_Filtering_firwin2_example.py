import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

fs, numtaps = 2, 41  # sampling frequency and number of taps
freqs = [0.0, 0.3, 0.6, 1.0]  # corner frequencies
gains = [1.0, 2.0, 0.5, 0.8]  # desired gains at corner frequencies

b = signal.firwin2(numtaps, freqs, gains, fs=fs)  # design filter
f, H = signal.freqz(b, fs=fs)  # calculate frequency response

fg, ax = plt.subplots(1, 1, layout="constrained")  # do the plotting
ax.set(title=f"Frequency Response of {numtaps}-tap Band-pass FIR-Filter",
       ylabel="Gain", xlim=(0, fs/2),
       xlabel=rf"Frequency $f\,$ in hertz (sampling frequency $f_S={fs}\,$Hz)")
ax.plot(freqs, gains, 'ko--', alpha=.5, label="Desired Gains")
ax.plot(f, np.abs(H), label="Response $|H(f)|$")
ax.grid()
ax.legend()
plt.show()
