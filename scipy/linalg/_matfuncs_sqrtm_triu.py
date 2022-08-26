# pythran export within_block_loop(float64[:,:], float64[:,:],
#                                  (int, int) list, intp)
# pythran export within_block_loop(complex128[:,:], complex128[:,:],
#                                  (int, int) list, intp)

def within_block_loop(R, T, start_stop_pairs, nblocks):
    for start, stop in start_stop_pairs:
        for j in range(start, stop):
            for i in range(j-1, start-1, -1):
                s = 0
                if j - i > 1:
                    # s = R[i, i+1:j] @ R[i+1:j, j]
                    for k in range(i + 1, j):
                        s += R[i, k] * R[k, j]

                denom = R[i, i] + R[j, j]
                num = T[i, j] - s
                if denom != 0:
                    R[i, j] = (T[i, j] - s) / denom
                elif denom == 0 and num == 0:
                    R[i, j] = 0
                else:
                    raise RuntimeError('failed to find the matrix square root')
