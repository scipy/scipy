import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multiscale_graphcorr


def mgc_plot(x, y, sim_name, mgc_dict):
    """Plot sim and MGC-plot"""
    # local correlation map
    plt.figure(figsize=(8, 8))
    mgc_map = mgc_dict["mgc_map"]

    ax = plt.gca()

    # draw heatmap
    plt.title("Local Correlation Map", fontsize=17)
    im = ax.imshow(mgc_map, cmap='YlGnBu')

    # colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("", rotation=-90, va="bottom")
    ax.invert_yaxis()

    # Turn spines off and create white grid.
    for _, spine in ax.spines.items():
        spine.set_visible(False)

    # optimal scale
    opt_scale = mgc_dict["opt_scale"]
    ax.scatter(opt_scale[0], opt_scale[1],
               marker='X', s=200, color='red')

    # other formatting
    ax.tick_params(bottom="off", left="off")
    ax.set_xlabel('#Neighbors for X', fontsize=15)
    ax.set_ylabel('#Neighbors for Y', fontsize=15)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    plt.show()


x = np.array([-0.57765724, -2.9276214, -1.04262652, 1.2582997, 2.02245107,
              0.17284193, 0.32260546, 3.57466422, -4.48283493, -1.54766131,
              1.04438668, -2.75866533, -0.69001856, -2.88075929, 3.64557631,
              2.5336397, 3.25705372, -0.45880268, 0.14090932, 1.99165555,
              -1.14072871, 3.0244705, 3.99267301, 0.58097879, 1.91957596,
              -0.09198175, -0.0202113, -0.53583113, 1.04351052, -0.94095353,
              1.25420529, 1.82505601, 1.76087592, -0.84837026, -0.35739805,
              -0.55455663, -0.98819013, 0.23788301, -1.00355776, 3.34646196,
              2.82813969, 3.75906558, -1.35422137, -0.86101343, -0.78098675,
              1.29833785, -0.51161473, -0.6701878, 2.94804317, -1.0313851,
              -2.55195402, 0.17734692, -0.0494121, 0.24185469, 1.83039907,
              3.75957814, 1.18814491, -2.25471647, -0.84504751, -0.67927423])


y = np.array([0.51031563, 0.38642023, -0.51907824, -1.27170782, 0.14132488,
              0.20053897, 2.38884352, -1.49329281, 1.88817141, 4.16704358,
              -1.46519457, -1.53420917, -1.05129013, 0.64661363, -1.20961109,
              3.63823477, 2.67217981, 0.54437479, 0.29491934, 0.68780628,
              -2.8710651, 2.86781294, 0.62258443, -1.74076219, -0.18287893,
              -1.61226611, 0.58488943, 0.38614018, 2.03597632, -0.73118874,
              -3.45631826, -0.51425533, -3.13829749, -1.01595659, -3.39980512,
              -3.46500481, -0.23870139, -1.49681269, -0.71348741, -1.88981505,
              -2.61335135, -0.90048746, 2.3437962, -1.11187083, -3.42930925,
              -1.28428743, 0.60948834, 0.47984385, 3.10223942, -0.12788552,
              1.32301543, 0.23418854, 0.49738729, -1.52406425, -0.5740983,
              1.72467609, -1.32066007, 1.84240793, -1.00591703, 0.29352814])


_, _, mgc_dict = multiscale_graphcorr(x, y)
mgc_plot(x, y, "Spiral", mgc_dict)
