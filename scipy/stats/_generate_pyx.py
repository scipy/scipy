import pathlib

def make_biasedurn():
    # numpy's random C API changes between 17 and 18
    from numpy import __version__ as __numpy_version__
    from scipy._lib import _pep440
    NPY_OLD = _pep440.parse(__numpy_version__) < _pep440.Version('1.18')
    del _pep440
    biasedurn_base = (pathlib.Path(__file__).parent / 'biasedurn').absolute()
    with open(biasedurn_base.with_suffix('.pyx.templ'), 'r') as src:
        contents = src.read()
    with open(biasedurn_base.with_suffix('.pyx'), 'w') as dest:
        dest.write(contents.format(NPY_OLD=str(bool(NPY_OLD))))


if __name__ == '__main__':
    make_biasedurn()
