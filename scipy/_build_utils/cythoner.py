#!python3
import os
import os.path as op
import sys
import subprocess as sbp

from Cython.Compiler.Main import Context

_, in_fname, out_fname = sys.argv
in_fname, out_fname = (op.abspath(p) for p in (in_fname, out_fname))

module_name = Context([], []).extract_module_name(in_fname, object())

res = sbp.run(
    ['cython', '-3',
    '--module-name', module_name,
    '--output-file', out_fname,
    '--include-dir', os.getcwd(),
     in_fname],
    check=True)
