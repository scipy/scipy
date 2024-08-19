import matplotlib.pyplot as plt
import numpy as np

aa = [1, 1, 2, 2]  # amplitudes
taus = [1, 2, 1, 2]  # durations

fg0, axx = plt.subplots(3, 1, sharex='all', tight_layout=True, figsize=(6., 4.))
axx[0].set(title=r"Spectrum $|X(f)|$", ylabel="V/Hz")
axx[1].set(title=r"Magnitude Spectrum $|X(f)/\tau|$ ", ylabel=r"V")
axx[2].set(title=r"Amplitude Spectral Density $|X(f)/\sqrt{\tau}|$",
           ylabel=r"$\operatorname{V} / \sqrt{\operatorname{Hz}}$",
           xlabel="Frequency $f$ in Hertz",)

x_labels, x_ticks = [], []
f = np.linspace(-2.5, 2.5, 400)
for c_, (a_, tau_) in enumerate(zip(aa, taus), start=1):
    aZ_, f_ = abs(a_ * tau_ * np.sinc(tau_ * f) / 2), f + c_ * 5
    axx[0].plot(f_, aZ_)
    axx[1].plot(f_, aZ_ / tau_)
    axx[2].plot(f_, aZ_ / np.sqrt(tau_))
    x_labels.append(rf"$a={a_:g}$, $\tau={tau_:g}$")
    x_ticks.append(c_ * 5)

axx[2].set_xticks(x_ticks)
axx[2].set_xticklabels(x_labels)
plt.show()
