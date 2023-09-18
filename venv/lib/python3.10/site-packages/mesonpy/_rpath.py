# SPDX-FileCopyrightText: 2023 The meson-python developers
#
# SPDX-License-Identifier: MIT

from __future__ import annotations

import os
import subprocess
import sys
import typing


if typing.TYPE_CHECKING:
    from typing import List

    from mesonpy._compat import Iterable, Path


if sys.platform == 'linux':

    def _get_rpath(filepath: Path) -> List[str]:
        r = subprocess.run(['patchelf', '--print-rpath', os.fspath(filepath)], capture_output=True, text=True)
        return r.stdout.strip().split(':')

    def _set_rpath(filepath: Path, rpath: Iterable[str]) -> None:
        subprocess.run(['patchelf','--set-rpath', ':'.join(rpath), os.fspath(filepath)], check=True)

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        rpath = _get_rpath(filepath)
        if '$ORIGIN/' in rpath:
            rpath = [('$ORIGIN/' + libs_relative_path if path == '$ORIGIN/' else path) for path in rpath]
            _set_rpath(filepath, rpath)


elif sys.platform == 'darwin':

    def _get_rpath(filepath: Path) -> List[str]:
        rpath = []
        r = subprocess.run(['otool', '-l', os.fspath(filepath)], capture_output=True, text=True)
        rpath_tag = False
        for line in [x.split() for x in r.stdout.split('\n')]:
            if line == ['cmd', 'LC_RPATH']:
                rpath_tag = True
            elif len(line) >= 2 and line[0] == 'path' and rpath_tag:
                rpath.append(line[1])
                rpath_tag = False
        return rpath

    def _replace_rpath(filepath: Path, old: str, new: str) -> None:
        subprocess.run(['install_name_tool', '-rpath', old, new, os.fspath(filepath)], check=True)

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        rpath = _get_rpath(filepath)
        if '@loader_path/' in rpath:
            _replace_rpath(filepath, '@loader_path/', '@loader_path/' + libs_relative_path)

else:

    def fix_rpath(filepath: Path, libs_relative_path: str) -> None:
        raise NotImplementedError(f'Bundling libraries in wheel is not supported on {sys.platform}')
