#!python3
import os
import os.path as op
import sys
import subprocess as sbp

_, in_fname, out_fname = sys.argv
in_fname, out_fname = (op.abspath(p) for p in (in_fname, out_fname))

res = sbp.run(
    ['cython', '-3', '--fast-fail',
    '--output-file', out_fname,
    '--include-dir', os.getcwd(),
     in_fname],
    check=True)
