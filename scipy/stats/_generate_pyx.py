import pathlib

def make_biasedurn():
    # numpy's random C API changes between 1.17 and 1.18
    import numpy as np
    ver = np.__version__.split('.')
    NPY_OLD = (int(ver[0]) <= 1) and (int(ver[1]) < 18)
    biasedurn_base = (pathlib.Path(__file__).parent / 'biasedurn').absolute()
    with open(biasedurn_base.with_suffix('.pyx.templ'), 'r') as src:
        contents = src.read()
    with open(biasedurn_base.with_suffix('.pyx'), 'w') as dest:
        dest.write(contents.format(NPY_OLD=str(bool(NPY_OLD))))


if __name__ == '__main__':
    make_biasedurn()
