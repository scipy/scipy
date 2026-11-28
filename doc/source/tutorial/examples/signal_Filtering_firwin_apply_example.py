import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

fs = 100  # sampling frequency in Hz
t = np.arange(400) / fs  # time axis of the signal

# A low-frequency component of interest buried under higher-frequency tones:
clean = np.cos(2 * np.pi * 0.5 * t)  # 0.5 Hz signal to recover
x = (clean + 0.4 * np.sin(2 * np.pi * 15.3 * t)
     + 0.3 * np.sin(2 * np.pi * 23.45 * t + 0.8))

# Design a low-pass FIR with a Kaiser window. ``kaiserord`` returns the number
# of taps and the ``beta`` parameter needed to reach a 60 dB stop-band
# attenuation with a 5 Hz wide transition band:
nyq = fs / 2  # Nyquist frequency
numtaps, beta = signal.kaiserord(ripple=60, width=5 / nyq)
f_c = 3  # cut-off frequency (-6 dB gain) in Hz
taps = signal.firwin(numtaps, f_c, window=('kaiser', beta), fs=fs)

# Apply the filter. A linear-phase FIR filter of ``numtaps`` taps delays the
# signal by a constant group delay of ``(numtaps - 1) / 2`` samples, which is
# removed here so the filtered signal lines up with the input:
y = signal.lfilter(taps, 1.0, x)
delay = 0.5 * (numtaps - 1) / fs  # group delay in seconds

# do the plotting:
fg, ax = plt.subplots(1, 1, layout="constrained", figsize=(6, 4))
ax.set(title=f"Low-pass Filtering with a {numtaps}-tap FIR Filter",
       xlabel=rf"Time $t\,$ in seconds (sampling frequency $f_S={fs}\,$Hz)",
       ylabel="Amplitude", xlim=(t[0], t[-1]))
ax.plot(t, x, 'C0', alpha=0.4, label="Noisy signal $x(t)$")
ax.plot(t - delay, y, 'C3', label="Filtered signal (delay compensated)")
ax.plot(t, clean, 'C2--', alpha=0.7, label="0.5 Hz component")
ax.legend()
ax.grid()
plt.show()
