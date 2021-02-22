'''Helper functions to get location of header files.'''

import pathlib
from typing import List, Union


def _boost_dir(ret_path: bool = False) -> Union[pathlib.Path, str]:
    '''Directory where root Boost/ directory lives.'''
    p = pathlib.Path(__file__).parent / 'boost'
    return p if ret_path else str(p)


def _boost_inc_dirs(ret_path: bool = False) -> Union[List[pathlib.Path],
                                                     List[str]]:
    '''Root directories for all Boost headers.'''
    inc_dirs = []
    for lib in (_boost_dir(ret_path=True) / 'libs').glob('*'):
        if lib.is_dir():
            if lib.name == 'headers':
                continue
            if lib.name == 'numeric':
                for sublib in lib.glob('*'):
                    if sublib.is_dir() and (sublib / 'include').exists():
                        p = sublib / 'include'
                        inc_dirs.append(p if ret_path else str(p))
            else:
                p = lib / 'include'
                inc_dirs.append(p if ret_path else str(p))
    return inc_dirs
