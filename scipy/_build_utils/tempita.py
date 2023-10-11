import sys
import os
import argparse

from Cython import Tempita as tempita
# XXX: If this import ever fails (does it really?), vendor either
# cython.tempita or numpy/npy_tempita.


def process_tempita(fromfile, outfile=None):
    """Process tempita templated file and write out the result.

    The template file is expected to end in `.c.in` or `.pyx.in`:
    E.g. processing `template.c.in` generates `template.c`.

    """
    from_filename = tempita.Template.from_filename
    template = from_filename(fromfile,
                             encoding=sys.getdefaultencoding())

    content = template.substitute()

    with open(outfile, 'w') as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Path to the input file")
    parser.add_argument("-o", "--outdir", type=str,
                        help="Path to the output directory")
    parser.add_argument("-i", "--ignore", type=str,
                        help="An ignored input - may be useful to add a "
                             "dependency between custom targets")
    args = parser.parse_args()

    if not args.infile.endswith('.in'):
        raise ValueError(f"Unexpected extension: {args.infile}")

    if not args.outdir:
        raise ValueError("Missing `--outdir` argument to tempita.py")

    outdir_abs = os.path.join(os.getcwd(), args.outdir)
    outfile = os.path.join(outdir_abs,
                           os.path.splitext(os.path.split(args.infile)[1])[0])

    process_tempita(args.infile, outfile)


if __name__ == "__main__":
    main()
