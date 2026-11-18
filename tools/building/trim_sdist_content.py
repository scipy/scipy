#!/usr/bin/env python3
"""
The purpose of this script is to remove files from the sdist that are not
needed and bloat the sdist size too much. This deals with files from
git submodules, because those cannot be removed by using `export-ignore`
in the top-level `.gitattributes` file.
"""

import os
import pathlib
import shutil

dist_root = pathlib.Path(os.environ['MESON_DIST_ROOT'])

for name in [dist_root / d for d in (
    'subprojects/boost_math/math/.github',
    'subprojects/boost_math/math/build',
    'subprojects/boost_math/math/config',
    'subprojects/boost_math/math/doc',
    'subprojects/boost_math/math/example',
    'subprojects/boost_math/math/meta',
    'subprojects/boost_math/math/reporting',
    'subprojects/boost_math/math/src',
    'subprojects/boost_math/math/test',
    'subprojects/boost_math/math/tools',
    'subprojects/highs/.github',
    'subprojects/highs/app',
    'subprojects/highs/check',
    'subprojects/highs/docs',
    'subprojects/highs/examples',
    'subprojects/highs/nuget',
    'subprojects/highs/scripts',
    'subprojects/highs/tests',
    'subprojects/xsf/.github',
    'subprojects/xsf/pixi.lock',
    'subprojects/xsf/tests',
    )]:
    if name.is_file():
        name.unlink()
    else:
        shutil.rmtree(name)

