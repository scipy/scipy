#!/usr/bin/env python

import re
from itertools import chain
from glob import iglob
import sys
import argparse
import os
from get_submodule_paths import get_submodule_paths


# The set of Unicode code points greater than 127 that we allow in the source code:
# Note that the lines enclosed by the marker lines BEGIN_INCLUDE_RST/END_INCLUDE_RST
# are included in the file `doc/source/dev/missing-bits.rst` to be rendered in the user
# documentation.
#
# BEGIN_INCLUDE_RST (do not change this line!)
latin1_letters = set(chr(cp) for cp in range(192, 256))
greek_letters = set('αβγδεζηθικλμνξoπρστυϕχψω' + 'ΓΔΘΛΞΠΣϒΦΨΩ')
box_drawing_chars = set(chr(cp) for cp in range(0x2500, 0x2580))
extra_symbols = set('®ő∫≠≥≤±∞²³·→√✅⛔⚠️')
allowed = latin1_letters | greek_letters | box_drawing_chars | extra_symbols
# END_INCLUDE_RST (do not change this line!)


def check_unicode(showall=False):
    """
    If showall is True, all non-ASCII characters are displayed.
    """
    # File encoding regular expression from PEP-263.
    encoding_pat = re.compile("^[ \t\f]*#.*?coding[:=][ \t]*([-_.a-zA-Z0-9]+)")
    root_dir = os.path.dirname(os.path.dirname(__file__))
    submodule_paths = get_submodule_paths()

    nbad = 0
    for name in chain(iglob(os.path.join(root_dir, 'scipy/**/*.py'), recursive=True),
                      iglob(os.path.join(root_dir, 'scipy/**/*.pyx'), recursive=True),
                      iglob(os.path.join(root_dir, 'scipy/**/*.px[di]'),
                            recursive=True)):
        if any(submodule_path in name for submodule_path in submodule_paths):
            continue
        # Read the file as bytes, and check for any bytes greater than 127.
        with open(name, 'rb') as f:
            content = f.read()
        if len(content) == 0:
            continue
        if max(content) > 127:
            # There is at least one non-ASCII character in the file.
            # Check the first two lines for an encoding comment.
            lines = content.splitlines()
            for line in lines[:2]:
                match = re.match(encoding_pat,
                                 line.decode(encoding='latin-1'))
                if match:
                    break

            # If an explicit encoding was given in a comment, use
            # that to decode the contents. Otherwise use UTF-8.
            if match:
                encoding = match[1]
                file_enc_msg = f"(explicit encoding '{encoding}')"
            else:
                encoding = 'utf-8'
                file_enc_msg = "(no explicit encoding; utf-8 assumed)"
            content = content.decode(encoding=encoding)

            out = []
            for n, line in enumerate(content.splitlines()):
                for pos, char in enumerate(line):
                    cp = ord(char)
                    if cp > 127:
                        msg = (f"... line {n+1}, position {pos+1}: "
                               f"character '{char}', code point U+{cp:04X}")
                        if showall:
                            out.append(msg)
                        else:
                            if char not in allowed:
                                out.append(msg)
            if len(out) > 0:
                nbad += 1
                print(f"{name} {file_enc_msg}")
                for msg in out:
                    print(msg)
    return nbad


if __name__ == "__main__":
    descr = ('Check for disallowed Unicode characters in the SciPy Python and '
             ' Cython source code.')
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--showall', action='store_true',
                        help=('Show non-ASCII Unicode characters from all '
                              'files.'))
    args = parser.parse_args()
    sys.exit(check_unicode(args.showall) > 0)
