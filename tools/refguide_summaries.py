#!/usr/bin/env python
"""Generate function summaries for the refguide. For example, if the
__init__ file of a submodule contains:

.. autosummary::
   :toctree: generated/

   foo
   foobar

Then it will modify the __init__ file to contain (*)

.. autosummary::
   :toctree: generated/

   foo    -- First line of the documentation of `foo`.
   foobar -- First line of the documentation of `foobar`.

If there is already text after the function definitions it will be
overwritten, i.e.

.. autosummary::
   :toctree: generated/

   foo    -- Blah blah blah.
   foobar -- Blabbity blabbity.

will also become (*).

"""
import os
import argparse
import importlib
import re


EXCEPTIONS = {
    'jn': ('Bessel function of the first kind of real order and '
           'complex argument')
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("module",
                        help="module to add summaries to")
    parser.add_argument("--dry-run",
                        help="print __init__ file instead of overwriting",
                        action="store_true")
    args = parser.parse_args()

    filename = os.path.join(os.path.dirname(__file__), '..', 'scipy',
                            args.module, '__init__.py')
    module = importlib.import_module('scipy.' + args.module)

    fnew = []
    with open(filename) as f:
        line = f.readline()
        while line:
            if '.. autosummary::' in line:
                fnew.append(line.rstrip())
                fnew.append(f.readline().rstrip())  # :toctree: generated/
                fnew.append(f.readline().rstrip())  # blank line
                line = f.readline()
                summaries = []
                maxlen = 0
                while line.strip():
                    func = line.split('--')[0].strip()
                    ufunc = '[+]' not in line
                    if len(func) > maxlen:
                        maxlen = len(func)

                    if func in EXCEPTIONS.keys():
                        summary = [EXCEPTIONS[func]]
                    else:
                        summary = []
                        doc = getattr(module, func).__doc__.split('\n')
                        i = 0 if doc[0].strip() else 1
                        while True:
                            if re.match(func + r'\(.*\)', doc[i].strip()):
                                # ufunc docstrings contain the signature
                                i += 2
                            else:
                                break
                        while i < len(doc) and doc[i].strip():
                            summary.append(doc[i].lstrip())
                            i += 1

                    summary = ' '.join([x.lstrip() for x in summary])
                    summary = '[+]' + summary if not ufunc else summary
                    summaries.append((func, summary))
                    line = f.readline()
                for (func, summary) in summaries:
                    spaces = ' '*(maxlen - len(func) + 1)
                    fnew.append('   ' + func + spaces + '-- ' + summary)
                fnew.append(line.rstrip())
            else:
                fnew.append(line.rstrip())
            line = f.readline()

    if args.dry_run:
        print('\n'.join(fnew))
    else:
        with open(filename, 'w') as f:
            f.write('\n'.join(fnew))
            f.write('\n')


if __name__ == "__main__":
    main()
