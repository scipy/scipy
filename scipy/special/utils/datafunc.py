import csv

import numpy as np


def parse_txt_data(filename):
    f = open(filename)
    try:
        reader = csv.reader(f, delimiter=',')
        data = [list(map(float, row)) for row in reader]
        nc = len(data[0])
        for i in data:
            if not nc == len(i):
                raise ValueError(i)
        ## guess number of columns/rows
        #row0 = f.readline()
        #nc = len(row0.split(',')) - 1
        #nlines = len(f.readlines()) + 1
        #f.seek(0)
        #data = np.fromfile(f, sep=',')
        #if not data.size == nc * nlines:
        #    raise ValueError("Inconsistency between array (%d items) and "
        #                     "guessed data size %dx%d" % (data.size, nlines, nc))
        #data = data.reshape((nlines, nc))
        #return data
    finally:
        f.close()

    return np.array(data)


def run_test(filename, funcs, args=[0]):  # noqa: B006
    nargs = len(args)
    if len(funcs) > 1 and nargs > 1:
        raise ValueError("nargs > 1 and len(funcs) > 1 not supported")

    data = parse_txt_data(filename)
    if data.shape[1] != len(funcs) + nargs:
        raise ValueError("data has %d items / row, but len(funcs) = %d and "
                         "nargs = %d" % (data.shape[1], len(funcs), nargs))

    if nargs > 1:
        f = funcs[0]
        x = [data[args[i]] for i in nargs]
        return f(*x)
    else:
        y = [f(data[:, 0]) - data[:, idx + 1] for idx, f in enumerate(funcs)]

        return data[:, 0], y


if __name__ == '__main__':
    from convert import DATA_DIR
    import os

    data = []
    for root, dirs, files in os.walk(DATA_DIR):
        for f in files:
            name = os.path.join(root, f)
            print(name)
            data.append(parse_txt_data(name))
