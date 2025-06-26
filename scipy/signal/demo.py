# Demo file for plotting modulated amplitude wave

import numpy as np
import matplotlib.pyplot as plt
t = np.linspace(0,1,1000)

carrier = np.sin(2*np.pi*50*t) # 50 Hz wave
modulator = 1 + 0.5*np.sin(2*np.pi*5*t) # 5 Hz wave
modulated_signal = carrier * modulator

#Plotting
plt.plot(t, modulated_signal)
plt.title('Amplitude Modulated Signal')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.grid(True)
plt.show()