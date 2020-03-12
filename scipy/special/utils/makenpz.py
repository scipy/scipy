"""
python makenpz.py DIRECTORY

Build a npz containing all data files in the directory.

"""

import os
import numpy as np
import argparse

from distutils.util import newer


def main():
    p = argparse.ArgumentParser(usage=(__doc__ or '').strip())
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
    np.savez_compressed(outp, **data)


if __name__ == "__main__":
    main()
