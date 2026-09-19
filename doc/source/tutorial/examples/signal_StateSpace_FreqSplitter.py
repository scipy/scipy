import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

ABCD = (np.array([[-10*np.sqrt(2), -100], [1,   0]]), np.array([[1], [0]]),
        np.array([[-10*np.sqrt(2), -100], [0, 100]]), np.array([[1], [0]]))

t_i, yy_i = signal.impulse(ABCD)  # simulate impulse response
t_s, yy_s = signal.step(ABCD)    # simulate step response

(b0, b1), a = signal.ss2tf(*ABCD)  # convert to transfer functions
f = np.geomspace(1e-1, 1e2, 100)  # log-spaced frequencies

_, (ax0, ax1) = plt.subplots(2, 1, sharex='all', constrained_layout=True)
ax0.set(title="Impulse response", ylabel="Amplitude", xlim=(t_i[0], t_i[-1]))
ax1.set(title="Step response", ylabel="Amplitude", xlabel="Time in seconds")
for c_ in range(2):
    ax0.plot(t_i, yy_i[:, c_], f'C{c_}-', label=f"Output {c_}")
    ax1.plot(t_s, yy_s[:, c_], f'C{c_}-', label=f"Output {c_}")

_, (ax10, ax11) = plt.subplots(2, 1, sharex='all', constrained_layout=True)
ax10.set_title("Frequency response")
ax10.set(ylabel="Magnitude in dB", xscale='log', xlim=(f[0], f[-1]),
         yticks=np.arange(-80, 20, 20))
ax11.set(ylabel="Phase in radians", xlabel="Frequency in hertz",
         yticks = np.pi*np.array([-1, -0.5, 0, 0.5, 1]),
         yticklabels=['-π', '-π/2', '0', 'π/2', 'π'])
for c_, b_ in enumerate((b0, b1)):
    _, H = signal.freqs(b_, a, worN=f)  # calculate frequency response
    H_dB, H_ang = 20*np.log10(np.abs(H)), np.unwrap(np.angle(H))
    ax10.plot(f, H_dB, f'C{c_}-', label=f"Output {c_}")
    ax11.plot(f, H_ang, f'C{c_}-', label=f"Output {c_}")

for ax_ in (ax0, ax1, ax10, ax11):
    ax_.grid(True)
    ax_.legend()
plt.show()



