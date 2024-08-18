import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import rfft, rfftfreq

n, T = 100, 0.01  # number of samples and sampling interval
tau = n*T
t = np.arange(n) * T
fcc = (20, 20.5)  # frequencies of sines
xx = (np.sin(2 * np.pi * fc_ * t) for fc_ in fcc)  # sine signals

f = rfftfreq(n, T)  # frequency bins range from 0 Hz to Nyquist freq.
XX = (rfft(x_) / n for x_ in xx)  # one-sided FFT normalized to magnitude

i0, i1 = 15, 25
f_cont = np.linspace(f[i0], f[i1], 501)

fg1, axx = plt.subplots(1, 2, sharey='all', tight_layout=True,
                        figsize=(6., 3.))
for c_, (ax_, X_, fx_) in enumerate(zip(axx, XX, fcc)):
    Xc_ = (np.sinc(tau * (f_cont - fx_)) +
           np.sinc(tau * (f_cont + fx_))) / 2
    ax_.plot(f_cont, abs(Xc_), f'-C{c_}', alpha=.5, label=rf"$f_x={fx_}\,$Hz")
    m_line, _, _, = ax_.stem(f[i0:i1+1], abs(X_[i0:i1+1]), markerfmt=f'dC{c_}',
                             linefmt=f'-C{c_}', basefmt=' ')
    plt.setp(m_line, markersize=5)

    ax_.legend(loc='upper left', frameon=False)
    ax_.set(xlabel="Frequency $f$ in Hertz", xlim=(f[i0], f[i1]),
            ylim=(0, 0.59))

axx[0].set(ylabel=r'Magnitude $|X(f)/\tau|$')
fg1.suptitle("Continuous and Sampled Magnitude Spectrum ", x=0.55, y=0.93)
fg1.tight_layout()
plt.show()
