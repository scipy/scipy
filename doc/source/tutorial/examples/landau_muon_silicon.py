from functools import partial
from scipy import stats
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import numpy as np


def landau_energy_loss(E, E_mpv, xi):
    """Landau energy loss PDF

    Args:
        E: Energy loss variable
        E_mpv: Most probable energy loss
        xi: Width parameter
    """
    landau_loc = E_mpv + xi * (1 - np.euler_gamma - 0.20005183774398613)
    return stats.landau.pdf(
        E, loc=landau_loc + xi * np.log(np.pi / 2), scale=xi * np.pi / 2
    )


# Most probable energy loss for 10 GeV muon in 1.7 mmm of Silicon
# See PDG 2025, Fig. 34.7
E_mpv = 0.5268
# Width parameter (~1/4 of the FWHM)
# See PDG 2025, following Eqn. 34.12
xi = 0.03031
loss = partial(landau_energy_loss, E_mpv=E_mpv, xi=xi)

evals = np.linspace(0.4, 1.0, 400)
dist_scipy = loss(evals)

lmax = loss(E_mpv)
lo_halfmax = root_scalar(lambda E: loss(E) - lmax / 2, bracket=(0.4, E_mpv))
hi_halfmax = root_scalar(lambda E: loss(E) - lmax / 2, bracket=(E_mpv, 1.0))
fwhm = hi_halfmax.root - lo_halfmax.root


fig, ax = plt.subplots(1, 1)
ax.plot(evals, dist_scipy, label="Landau")
ax.axvline(E_mpv, ls=":", label="Mode")
ax.plot(
    [lo_halfmax.root, hi_halfmax.root],
    [lmax / 2, lmax / 2],
    label=rf"FWHM = ${fwhm / xi:.4f}\xi$",
)
ax.set_ylabel("Probability Density")
ax.set_xlabel("Energy Loss (MeV)")
ax.set_ylim(0, None)
ax.set_xlim(0.4, 1.0)
ax.set_xticks(np.arange(0.4, 1.0, 0.02), minor=True)
ax.legend(title="10 GeV muon in 1.7 mm Si") 
plt.show()
