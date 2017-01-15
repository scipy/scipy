#!/usr/bin/env bash
"""
makenpz.py DIRECTORY

Build a npz containing all data files in the directory.

"""

from __future__ import division, print_function, absolute_import

import os
import numpy as np
from optparse import OptionParser


def main():
    p = OptionParser()
    options, args = p.parse_args()

    if len(args) != 1:
        p.error("no valid directory given")

    inp = args[0]
    outp = inp + ".npz"

    files = []
    for dirpath, dirnames, filenames in os.walk(inp):
        for fn in filenames:
            if fn.endswith('.txt'):
                files.append(
                    (dirpath[len(inp)+1:] + '/' + fn[:-4],
                     os.path.join(dirpath, fn)))

    data = {}
    for key, fn in files:
        key = key.replace('/', '-').strip('-')
        try:
            data[key] = np.loadtxt(fn)
        except ValueError:
            print("Failed to load", fn)

    savez_compress(outp, **data)


def savez_compress(file, *args, **kwds):
    # Import is postponed to here since zipfile depends on gzip, an optional
    # component of the so-called standard library.
    import zipfile
    # Import deferred for startup time improvement
    import tempfile

    if isinstance(file, str):
        if not file.endswith('.npz'):
            file = file + '.npz'

    namedict = kwds
    for i, val in enumerate(args):
        key = 'arr_%d' % i
        if key in namedict:
            raise ValueError("Cannot use un-named variables and keyword %s" % key)
        namedict[key] = val

    zip = zipfile.ZipFile(file, mode="w", compression=zipfile.ZIP_DEFLATED)

    # Stage arrays in a temporary file on disk, before writing to zip.
    fd, tmpfile = tempfile.mkstemp(suffix='-numpy.npy')
    os.close(fd)
    try:
        for key, val in namedict.items():
            fname = key + '.npy'
            fid = open(tmpfile, 'wb')
            try:
                np.lib.format.write_array(fid, np.asanyarray(val))
                fid.close()
                fid = None
                zip.write(tmpfile, arcname=fname)
            finally:
                if fid:
                    fid.close()
    finally:
        os.remove(tmpfile)

    zip.close()


if __name__ == "__main__":
    main()
