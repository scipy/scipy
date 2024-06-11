import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import rfft, rfftfreq
from scipy.signal import get_window

n, n_zp = 128, 16384  # number of samples without and with zero-padding
t = np.arange(n)
f = rfftfreq(n_zp, 1 / n)

ww = ['boxcar', 'hann', 'hamming', 'tukey', 'blackman', 'flattop']
fg0, axx = plt.subplots(len(ww), 1, sharex='all', sharey='all', figsize=(6., 4.))
for c_, (w_name_, ax_) in enumerate(zip(ww, axx)):
    w_ = get_window(w_name_, n, fftbins=False)
    W_ = rfft(w_ / abs(sum(w_)), n=n_zp)
    W_dB = 20*np.log10(np.maximum(abs(W_), 1e-250))
    ax_.plot(f, W_dB, f'C{c_}-', label=w_name_)
    ax_.text(0.1, -50, w_name_, color=f'C{c_}', verticalalignment='bottom',
             horizontalalignment='left', bbox={'color': 'white', 'pad': 0})
    ax_.set_yticks([-20, -60])
    ax_.grid(axis='x')

axx[0].set_title("Spectral Leakage of various Windows")
fg0.supylabel(r"Normalized Magnitude $20\,\log_{10}|W(f)/c^\operatorname{amp}|$ in dB",
              x=0.04, y=0.5, fontsize='medium')
axx[-1].set(xlabel=r"Normalized frequency $f/\Delta f$ in bins",
            xlim=(0, 9), ylim=(-75, 3))

fg0.tight_layout(h_pad=0.4)
plt.show()
