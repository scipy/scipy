import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from scipy.fft import fft, fftshift

z_pts = [  0.1-27.8j,  12.5-48.9j,  40.3-60.8j,  69.0-36.8j,  77.1- 8.1j,  77.5+24.1j,
          62.2+54.2j,  37.8+67.4j,  41.7+56.2j,  53.1+29.6j,  55.3- 3.1j,  49.4-39.3j,
          24.2-47.3j,   7.8-23.3j,   9.7+ 8.8j,   5.6+40.9j,  -5.2+40.5j,  -9.5+ 8.5j,
          -7.7-23.5j, -24.4-47.7j, -49.6-39.2j, -55.6- 2.8j, -53.4+29.8j, -41.4+56.3j,
         -37.8+67.6j, -62.3+54.7j, -77.4+24.2j, -77.3- 7.9j, -69.1-36.4j, -40.2-60.9j,
         -12.4-48.9j]
z_pts = np.asarray(z_pts)  # the samples
cc = fftshift(fft(z_pts, norm='forward'))  # calculate Fourier coefficients

# Use Fourier series to calculate curve points:
N, k0 = len(z_pts), -(len(cc) // 2)  # number of samples, first Fourier series index
t = np.linspace(0, 1, N * 8, endpoint=False)
z = sum(c_k * np.exp(2j * np.pi * k * t) for k, c_k in enumerate(cc, start=k0))

fg0, ax0 = plt.subplots(1, 1, tight_layout=True)
ax0.set(title=f"Parametric Closed Curve $z(t)$ created from {N} Samples",
        xlabel=r"Re$\{z(t)\}$", ylabel=r"Im$\{z(t)\}$")

ax0.scatter(z_pts.real, z_pts.imag, c=np.arange(N) / N, s=8, cmap='twilight')

lpts = np.stack((z.real, z.imag), axis=1)  # extract 2d points
ln0 = ax0.add_collection(LineCollection(np.stack((lpts[:-1], lpts[1:]), axis=1),
                                        array=t[:-1], cmap='twilight'))
cb0 = fg0.colorbar(ln0, label="Parameter $t$")
plt.show()
