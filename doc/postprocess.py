#!/usr/bin/env python
"""
%prog MODE FILES...

Post-processes HTML and Latex files output by Sphinx.
MODE is either 'html' or 'tex'.

"""
import sys
import io
import re, optparse

def main():
    p = optparse.OptionParser(__doc__)
    options, args = p.parse_args()

    if len(args) < 1:
        p.error('no mode given')

    mode = args.pop(0)

    if mode not in ('html', 'tex'):
        p.error('unknown mode %s' % mode)

    def _open(fn, *args, **kwargs):
        r"""Handle UTF-8 encoding when loading under Py3"""
        # Issue is that default encoding under Py3 is
        # locale dependent (might even be ASCII),
        # so need to specify the encoding.  Py2 doesn't care.
        if sys.version_info.major < 3:
            return open(fn, *args, **kwargs)
        return io.open(fn, *args, encoding='utf-8', **kwargs)

    for fn in args:
        with _open(fn, 'r') as f:
            if mode == 'html':
                lines = process_html(fn, f.readlines())
            elif mode == 'tex':
                lines = process_tex(f.readlines())

        with _open(fn, 'w') as f:
            f.write("".join(lines))


def process_html(fn, lines):
    return lines


def process_tex(lines):
    """
    Fix autosummary LaTeX bug in Sphinx < 1.7.3
    (cf https://github.com/sphinx-doc/sphinx/issues/4790)

    """
    new_lines = []
    for line in lines:
        line = line.replace(r'p{0.5\linewidth}', r'\X{1}{2}')

        new_lines.append(line)
    return new_lines


if __name__ == "__main__":
    main()
