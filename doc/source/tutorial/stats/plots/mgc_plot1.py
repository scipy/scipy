import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multiscale_graphcorr


def sim_plot(x, y, sim_name):
    """Plot simulations."""
    fig = plt.figure(figsize=(8,8))
    fig.suptitle(sim_name + " Simulation", fontsize=17)
    plt.scatter(x, y)
    plt.show()


x = np.array([0.09762701, 0.43037873, 0.20552675, 0.08976637,
              -0.1526904, 0.29178823, -0.12482558, 0.783546,
              0.92732552, -0.23311696, 0.58345008, 0.05778984,
              0.13608912, 0.85119328, -0.85792788, -0.8257414,
              -0.95956321, 0.66523969, 0.5563135, 0.7400243])
y = np.array([-0.75550809, 1.40576643, -0.04929934, -0.12927078,
              -0.77908808, 0.6805334, -0.9317745, 0.67717586,
              0.47959224, -0.03966571, 0.32804751, -0.53252625,
              0.12199801, 1.06535921, -0.82466927, -0.67450545,
              -1.27672425, 0.48386911, 0.22008328, 0.56024772])


sim_plot(x, y, "Linear")
