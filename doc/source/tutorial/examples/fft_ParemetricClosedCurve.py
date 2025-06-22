import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from scipy.fft import ifft, fftshift

X = np.array([0, -3.699 - 3669j, -0.4812 - 2650j, 1.013 - 2008j, 1.98 + 1253j,
              3.538 + 505.4j, -0.742 + 616.5j, 0, 2.91 + 236.9j, 1.498 + 187j,
              -4.823 + 305j, 0, 0, 0, 0, 0, 0, 0, 0, 5.863 + 165.6j, 0, 0,
              -3.965 + 146.9j, 0, 0.5565 + 290.3j, 0, -2.706 - 615.6j, 5.784 + 1870j,
              -6.915 - 1080j, 8.614 - 864j, 4.133 + 2527j]) # DFT values
N, N1 = len(X), len(X)*8  # number of samples and number of interpolation points

x = ifft(X)  # determine sample values of closed curve

t1 = np.arange(N1 + 1) / N1  # interpolate by utilizing Fourier series
x1 = sum([c_k * np.exp(2j * np.pi * k_ * t1)
          for k_, c_k in enumerate(fftshift(X / N), start=-(N//2))])

fg0, ax0 = plt.subplots(1, 1, tight_layout=True)
ax0.set(title=f"Parametric Closed Curve $x(t)$ generated from {len(x)} Samples",
        xlabel=r"Re$\{x(t)\}$", ylabel=r"Im$\{x(t)\}$")

ax0.scatter(x.real, x.imag, c=np.arange(N) / N, s=8, cmap='twilight')

pts1 = np.stack((x1.real, x1.imag), axis=1)  # extract 2d points
ln1 = ax0.add_collection(LineCollection(np.stack((pts1[:-1], pts1[1:]), axis=1),
                                        array=t1[:-1], cmap='twilight'))
fg0.colorbar(ln1, label="Parameter $t$")

plt.show()
