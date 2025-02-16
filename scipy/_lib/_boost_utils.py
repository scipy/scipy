'''Helper functions to get location of header files.'''

import pathlib


def _boost_dir(ret_path: bool = False) -> pathlib.Path | str:
    '''Directory where root Boost/ directory lives.'''
    p = pathlib.Path(__file__).parent / 'boost_math/include'
    return p if ret_path else str(p)
