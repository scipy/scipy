#pythran export _sosfilt_float(float[:, :], float[:, :], float[:, :, :])
def _sosfilt_float(sos, x, zi):
    # Modifies x and zi in place
    n_signals = x.shape[0]
    n_samples = x.shape[1]
    n_sections = sos.shape[0]

    for i in range(n_signals):
        zi_slice = zi[i, :, :]
        for n in range(n_samples):

            x_cur = x[i, n]

            for s in range(n_sections):
                x_new = sos[s, 0] * x_cur + zi_slice[s, 0]
                zi_slice[s, 0] = (sos[s, 1] * x_cur - sos[s, 4] * x_new
                                  + zi_slice[s, 1])
                zi_slice[s, 1] = sos[s, 2] * x_cur - sos[s, 5] * x_new
                x_cur = x_new

            x[i, n] = x_cur
    return x, zi

