#!/usr/bin/env python
import os, subprocess

def main():
    os.chdir('sparsetools')
    files = ['csr.i', 'csc.i', 'coo.i', 'dia.i', 'bsr.i', 'csgraph.i']
    for fn in files:
        ret = subprocess.call(['swig', '-c++', '-python', fn])
        if ret != 0:
            raise RuntimeError("swig failed")

if __name__ == "__main__":
    main()
