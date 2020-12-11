import pathlib
from zipfile import ZipFile
from time import time

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('boostinator', parent_package, top_path)

    # unzip boost archive if not already done
    _curdir = pathlib.Path(__file__).parent
    if not (_curdir / 'boost').exists():
        print('[boost] boost directory not found, unzipping...')
        zipfile = _curdir / 'boost.zip'
        t0 = time()
        ZipFile(zipfile).extractall(path=_curdir / 'boost')
        print(f'[boost] Took {time() - t0} seconds to extract boost')
    else:
        print('[boost] boost directory found!')

    config.add_subpackage('boost_test_build')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
