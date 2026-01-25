import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

fs = 2  # sampling frequency
n0, f0_c = 40, 0.5  # number of taps and cut-off frequency for H_0(f)
n1, ff1_c = 41, [0.3, 0.8]  # number of taps and cut-off frequencies for H_1(f)

b0 = signal.firwin(n0, f0_c, fs=fs)  # design FIR filters
b1 = signal.firwin(n1, ff1_c, fs=fs)

f0, H0 = signal.freqz(b0, fs=fs)  # calculate frequency responses
f1, H1 = signal.freqz(b1, fs=fs)
H0_dB, H1_dB = (20 * np.log10(np.abs(H_)) for H_ in (H0, H1))  # convert to dB

# do the plotting:
fg0, (ax0, ax1) = plt.subplots(2, 1, sharex='all', layout="constrained",
                               figsize=(6, 4))
ax0.set(title=f"Frequency Response of {n0}-tap Low-pass FIR-Filter", xlim=(0, fs/2))
ax0.plot(f0, H0_dB, 'C0', label="Response $|H_0(f)|$")
ax0.axvline(f0_c, color='C2', linestyle='--', label=rf"$f_0={f0_c}\,$Hz")

ax1.set_title(f"Frequency Response of {n1}-tap Band-stop FIR-Filter")
ax1.set_xlabel(rf"Frequency $f\,$ in hertz (sampling frequency $f_S={fs}\,$Hz)")
ax1.plot(f1, H1_dB, 'C0', label="Response $|H_1(f)|$")
ax1.axvline(ff1_c[0], color='C2', linestyle='--', label=rf"$f_0={ff1_c[0]}\,$Hz")
ax1.axvline(ff1_c[1], color='C3', linestyle='--', label=rf"$f_1={ff1_c[1]}\,$Hz")

for ax_ in (ax0, ax1):
    ax_.set_ylabel("Gain in dB")
    ax_.legend()
    ax_.grid()
plt.show()
