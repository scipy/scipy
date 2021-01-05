import pathlib
import numpy as np

def make_biasedurn():
    # numpy's random C API changes between 17 and 18
    NPY_OLD = int(np.__version__.split('.')[1]) < 18
    biasedurn_base = (pathlib.Path(__file__).parent / 'biasedurn').absolute()
    with open(biasedurn_base.with_suffix('.pyx.templ'), 'r') as src:
        contents = src.read()
    with open(biasedurn_base.with_suffix('.pyx'), 'w') as dest:
        dest.write(contents.format(NPY_OLD=str(bool(NPY_OLD))))


if __name__ == '__main__':
    make_biasedurn()
