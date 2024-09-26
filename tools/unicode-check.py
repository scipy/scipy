#!/usr/bin/env python

import re
from itertools import chain
from glob import iglob
import sys
import argparse


# The set of Unicode code points greater than 127 that we
# allow in the source code.
latin1_letters = set(chr(cp) for cp in range(192, 256))
box_drawing_chars = set(chr(cp) for cp in range(0x2500, 0x2580))
extra_symbols = set(['®', 'ő', 'λ', 'π', 'ω', '∫', '≠', '≥', '≤', 'μ',
                     '±', '∞'])
allowed = latin1_letters | box_drawing_chars | extra_symbols


def unicode_check(showall=False):
    """
    If showall is True, all non-ASCII characters are displayed.
    """
    # File encoding regular expression from PEP-263.
    encoding_pat = re.compile("^[ \t\f]*#.*?coding[:=][ \t]*([-_.a-zA-Z0-9]+)")

    nbad = 0
    for name in chain(iglob('scipy/**/*.py', recursive=True),
                      iglob('scipy/**/*.pyx', recursive=True),
                      iglob('scipy/**/*.px[di]', recursive=True)):
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
    sys.exit(unicode_check(args.showall) > 0)
