"""Helper functions to get location of header files."""

import pathlib
from typing import Union


def _highs_dir(ret_path: bool = False) -> Union[pathlib.Path, str]:
    """Directory where root Boost/ directory lives."""
    p = pathlib.Path(__file__).parent / 'highs'
    return p if ret_path else str(p)
