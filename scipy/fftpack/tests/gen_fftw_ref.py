from subprocess import Popen, PIPE, STDOUT

import numpy as np

SZ = [2, 3, 4, 8, 12, 15, 16, 17, 32, 64, 128, 256, 512, 1024]

def gen_data(dt):
    arrays = {}

    if dt == np.double:
        pg = './fftw_double'
    elif dt == np.float32:
        pg = './fftw_single'
    else:
        raise ValueError("unknown: %s" % dt)
    # Generate test data using FFTW for reference
    for type in [1, 2, 3, 4]:
        arrays[type] = {}
        for sz in SZ:
            a = Popen([pg, str(type), str(sz)], stdout=PIPE, stderr=STDOUT)
            st = [i.strip() for i in a.stdout.readlines()]
            arrays[type][sz] = np.fromstring(",".join(st), sep=',', dtype=dt)

    return arrays

data = gen_data(np.float32)
filename = 'fftw_single_ref'

# Save ref data into npz format
d = {}
d['sizes'] = SZ
for type in [1, 2, 3, 4]:
    for sz in SZ:
        d['dct_%d_%d' % (type, sz)] = data[type][sz]
np.savez(filename, **d)
