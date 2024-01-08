import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import rfft, rfftfreq
from scipy.signal.windows import hann

n, T = 100, 0.01  # number of samples and sampling interval
tau = n*T
q = 3  # over-sampling factor
t = np.arange(n) * T
fcc = (20, 20.5)  # frequencies of sines
xx = [np.sin(2 * np.pi * fc_ * t) for fc_ in fcc]  # sine signals
w = hann(n)
c_w = abs(sum(w))  # normalize constant for window

f_X = rfftfreq(n, T)  # frequency bins range from 0 Hz to Nyquist freq.
XX = (rfft(x_ * w) / c_w for x_ in xx)  # one-sided amplitude spectrum
# Oversampled spectrum:
f_Y = rfftfreq(n*q, T)  # frequency bins range from 0 Hz to Nyquist freq.
YY = (rfft(x_ * w, n=q*n) / c_w for x_ in xx)  # one-sided magnitude spectrum

i0, i1 = 15, 25
j0, j1 = i0*q, i1*q

fg1, axx = plt.subplots(1, 2, sharey='all', tight_layout=True,
                        figsize=(6., 3.))
for c_, (ax_, X_, Y_, fx_) in enumerate(zip(axx, XX, YY, fcc)):
    ax_.plot(f_Y[j0:j1 + 1], abs(Y_[j0:j1 + 1]), f'.-C{c_}',
             label=rf"$f_x={fx_}\,$Hz")
    m_ln, s_ln, _, = ax_.stem(f_X[i0:i1 + 1], abs(X_[i0:i1 + 1]), basefmt=' ',
                              markerfmt=f'dC{c_}', linefmt=f'-C{c_}')
    plt.setp(m_ln, markersize=5)
    plt.setp(s_ln, alpha=0.5)

    ax_.legend(loc='upper left', frameon=False)
    ax_.set(xlabel="Frequency $f$ in Hertz", xlim=(f_X[15], f_X[25]),
            ylim=(0, 0.59))

axx[0].set(ylabel=r'Magnitude $|X(f)/\tau|$')
fg1.suptitle(r"Magnitude Spectrum (Hann window, $%d\times$oversampled)" % q,
             x=0.55, y=0.93)
plt.show()
