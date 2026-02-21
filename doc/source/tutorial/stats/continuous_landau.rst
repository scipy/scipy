
.. _continuous-landau:

Landau distribution
===================

A special case of Lévy-stable distributions with :math:`\alpha=1`
and :math:`\beta=1` and support :math:`-\infty < x < \infty`. The probability
density function is given by

.. math::

    f(x) = \frac{1}{\pi}\int_0^\infty \exp(-t \log t - xt)\sin(\pi t) dt

The differential entropy is 2.37263644000448182, and the moments are undefined.

Implementation: `scipy.stats.landau`

Energy loss of charged particles in matter
------------------------------------------

The physical process of kinetic energy loss in charged particles, as they pass
through matter and ionize electrons, was first characterized [1]_ by Lev Landau
using this distribution that bears his namesake. This distribution is of limited
use in practice because it is idealized in that: its support has no upper bound,
while physically the energy loss cannot exceed the total energy of the incident
particle, corrected by Vavilov [2]_; and that it assumes an idealized
probability of single scattering proportional to :math:`1/r^2`. The latter
idealization has had improvements made by numerous researchers, see e.g.
[3]_. For further discussion on the modeling of passage of particles through matter,
see [4]_, Section 34.2.9.

In typical use in the context of energy loss in matter, the Landau distribution
is parameterized by its most probable energy loss value ``E_mpv`` and a width
parameter ``xi`` (:math:`\Delta_p` and :math:`\xi`, respectively, in Ref. [4]_,
Eqn. 34.12) The ``loc`` and ``scale`` parameters of the implementation in
`scipy.stats.landau` are related to these parameters as demonstrated in the
``landau_energy_loss`` function of the example below, which predicts the
distribution of energy loss for a muon traveling through Silicon.

.. plot:: tutorial/examples/landau_muon_silicon.py


References
----------
.. [1] L. Landau, "On the energy loss of fast particles by ionization," J. Phys. (USSR) 8 (1944) 201-205.
    http://e-heritage.ru/Book/10093344
.. [2] P. Vavilov, "Ionization losses of high-energy heavy particles," Sov. Phys. JETP 5 (1957) 749-751.
    https://jetp.ras.ru/cgi-bin/e/index/e/5/4/p749?a=list
.. [3] H. Bichsel, "Straggling in thin silicon detectors," Rev. Mod. Phys. 60 (1988) 663-699.
    :doi:`10.1103/RevModPhys.60.663`
.. [4] S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024) and 2025 update.
    https://pdg.lbl.gov/2025/reviews/rpp2025-rev-passage-particles-matter.pdf