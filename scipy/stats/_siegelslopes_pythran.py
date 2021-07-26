import numpy as np


# pythran export siegelslopes(int[:] or float[:],
#                             float[:],
#                             str)
def siegelslopes(y, x, method):
    deltax = np.expand_dims(x, 1) - x
    deltay = np.expand_dims(y, 1) - y
    slopes, intercepts = [], []

    for j in range(len(x)):
        id_nonzero, = np.nonzero(deltax[j, :])
        slopes_j = np.empty(len(id_nonzero))
        for i in range(len(id_nonzero)):
            slopes_j[i] = deltay[j, id_nonzero[i]] / deltax[j, id_nonzero[i]]
        slopes_j = deltay[j, id_nonzero] / deltax[j, id_nonzero]
        medslope_j = np.median(slopes_j)
        slopes.append(medslope_j)
        if method == 'separate':
            z = y*x[j] - y[j]*x
            intercept_j = np.empty(len(id_nonzero))
            for i in range(len(id_nonzero)):
                intercept_j[i] = z[id_nonzero[i]] / deltax[j, id_nonzero[i]]
            medintercept_j = np.median(intercept_j)
            intercepts.append(medintercept_j)

    medslope = np.median(np.asarray(slopes))
    if method == "separate":
        medinter = np.median(np.asarray(intercepts))
    else:
        medinter = np.median(y - medslope*x)

    return medslope, medinter
