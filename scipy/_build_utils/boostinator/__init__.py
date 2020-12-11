'''Helper functions to get location of header files.'''

import pathlib


def get_include_dir() -> pathlib.Path:
    '''Directory where root boost/ directory lives.'''
    return pathlib.Path(__file__).parent
