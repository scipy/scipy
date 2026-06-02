# This is a python implementation of dfoxs.m,
# provided at https://github.com/POptUS/BenDFO
import numpy as np


def dfoxs(n, nprob, factor):
    x = np.zeros(n)

    if nprob == 1 or nprob == 2 or nprob == 3:  # Linear functions.
        x = np.ones(n)
    elif nprob == 4:  # Rosenbrock function.
        x[0] = -1.2
        x[1] = 1
    elif nprob == 5:  # Helical valley function.
        x[0] = -1
    elif nprob == 6:  # Powell singular function.
        x[0] = 3
        x[1] = -1
        x[2] = 0
        x[3] = 1
    elif nprob == 7:  # Freudenstein and Roth function.
        x[0] = 0.5
        x[1] = -2
    elif nprob == 8:  # Bard function.
        x[0] = 1
        x[1] = 1
        x[2] = 1
    elif nprob == 9:  # Kowalik and Osborne function.
        x[0] = 0.25
        x[1] = 0.39
        x[2] = 0.415
        x[3] = 0.39
    elif nprob == 10:  # Meyer function.
        x[0] = 0.02
        x[1] = 4000
        x[2] = 250
    elif nprob == 11:  # Watson function.
        x = 0.5 * np.ones(n)
    elif nprob == 12:  # Box 3-dimensional function.
        x[0] = 0
        x[1] = 10
        x[2] = 20
    elif nprob == 13:  # Jennrich and Sampson function.
        x[0] = 0.3
        x[1] = 0.4
    elif nprob == 14:  # Brown and Dennis function.
        x[0] = 25
        x[1] = 5
        x[2] = -5
        x[3] = -1
    elif nprob == 15:  # Chebyquad function.
        for k in range(n):
            x[k] = (k + 1) / (n + 1)
    elif nprob == 16:  # Brown almost-linear function.
        x = 0.5 * np.ones(n)
    elif nprob == 17:  # Osborne 1 function.
        x[0] = 0.5
        x[1] = 1.5
        x[2] = 1
        x[3] = 0.01
        x[4] = 0.02
    elif nprob == 18:  # Osborne 2 function.
        x[0] = 1.3
        x[1] = 0.65
        x[2] = 0.65
        x[3] = 0.7
        x[4] = 0.6
        x[5] = 3
        x[6] = 5
        x[7] = 7
        x[8] = 2
        x[9] = 4.5
        x[10] = 5.5
    elif nprob == 19:  # Bdqrtic.
        x = np.ones(n)
    elif nprob == 20:  # Cube.
        x = 0.5 * np.ones(n)
    elif nprob == 21:  # Mancino.
        for i in range(n):
            ss = 0
            for j in range(n):
                frac = (i + 1) / (j + 1)
                ss = ss + np.sqrt(frac) * (
                    (np.sin(np.log(np.sqrt(frac)))) ** 5
                    + (np.cos(np.log(np.sqrt(frac)))) ** 5
                )
            x[i] = -8.710996e-4 * ((i - 49) ** 3 + ss)
    elif nprob == 22:  # Heart8ls.
        x = np.asarray([-0.3, -0.39, 0.3, -0.344, -1.2, 2.69, 1.59, -1.5])
    else:
        print(f"unrecognized function number {nprob}")
        return None

    return factor * x
