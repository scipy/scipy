"""
python makenpz.py DIRECTORY

Build a npz containing all data files in the directory.

"""

from __future__ import division, print_function, absolute_import

import os
import numpy as np
import argparse

from distutils.util import newer


def main():
    p = argparse.ArgumentParser(usage=__doc__.strip())
    p.add_argument('--use-timestamp', action='store_true', default=False,
                   help="don't rewrite npz file if it is newer than sources")
    p.add_argument('dirname')
    args = p.parse_args()

    inp = os.path.normpath(args.dirname)
    outp = inp + ".npz"

    # Skip rebuilding if no sources
    if os.path.isfile(outp) and not os.path.isdir(inp):
        print("[makenpz] {} not rebuilt".format(outp))
        return

    # Find source files
    files = []
    for dirpath, dirnames, filenames in os.walk(inp):
        for fn in filenames:
            if fn.endswith('.txt'):
                key = dirpath[len(inp)+1:] + '-' + fn[:-4]
                key = key.strip('-')
                files.append((key, os.path.join(dirpath, fn)))

    # Check if changes required
    if args.use_timestamp and os.path.isfile(outp):
        try:
            old_data = np.load(outp)
            try:
                changed = set(old_data.keys()) != set(key for key, _ in files)
            finally:
                old_data.close()
        except (IOError, OSError):
            # corrupted file
            changed = True

        changed = changed or any(newer(fn, outp) for key, fn in files)
        changed = changed or newer(__file__, outp)
        if not changed:
            print("[makenpz] {} is already up to date".format(outp))
            return

    data = {}
    for key, fn in files:
        data[key] = np.loadtxt(fn)

    print("[makenpz] generating {}".format(outp))
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
