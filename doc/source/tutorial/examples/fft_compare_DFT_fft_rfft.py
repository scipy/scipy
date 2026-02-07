import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq, irfft, rfft, rfftfreq

# Create signals:
x0 = np.array([4, 1, 0, 1])
x1 = irfft(rfft(x0), n=5)  # Here: rfft(x1) == rfft(x0)
tau = 1  # signal duration in seconds

ll = np.arange(-3, 6)  # indices of DFT coefficients
_, ax01 = plt.subplots(3, 2, sharex='col', tight_layout=True, figsize=(6., 4.))
for x, axx in zip((x1, x0), ax01.T):
    N = len(x)  # number of samples
    T, delta_f = tau / N, 1 / tau  # sampling interval and frequency resolution
    kk = np.arange(0, N)  # sample indices of DFT

    X_DFT = np.array([sum(x * np.exp(-2j * np.pi * l_ * kk / N)) for l_ in ll])
    f_DFT = ll * delta_f  # frequencies of DFT

    XX = fft(x), fft(x), rfft(x)  # DFT values and FFT values coincide
    ff = np.arange(N) * delta_f, fftfreq(N, T), rfftfreq(N, T)

    t_strs = ["DFT", "fft / fftfreq", "rfft / rfftfreq"]  # do the plotting:
    for p, (ax, f, X, ts) in enumerate(zip(axx, ff, XX, t_strs)):
        lns = axx[p].stem(f_DFT, abs(X_DFT), linefmt='k-', markerfmt="k.", basefmt=' ')
        [ln_.set_alpha(0.1) for ln_ in lns]  # set transparency
        axx[p].stem(f, abs(X), linefmt=f'C{p}-', markerfmt=f'C{p}.', basefmt=' ')
        ax.set(title=f"{ts} for len(x) = {N}", yticks=[], ylim=(0, 6.5))
    axx[-1].set(xticks=f_DFT, xlim=(f_DFT[0]-0.5, f_DFT[-1]+0.5))
    axx[-1].set_xlabel(rf"Frequency $f$ in hertz ($\Delta f ={delta_f:g}\,$Hz)")
plt.show()



