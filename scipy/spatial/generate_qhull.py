#!/usr/bin/env python
import tempfile
import subprocess
import os
import sys
import re
import shutil

tmp_dir = tempfile.mkdtemp()
try:
    # Run Cython
    dst_fn = os.path.join(tmp_dir, 'qhull.c')
    ret = subprocess.call(['cython', '-o', dst_fn, 'qhull.pyx'])
    if ret != 0:
        sys.exit(ret)

    # Strip comments
    f = open(dst_fn, 'r')
    text = f.read()
    f.close()

    r = re.compile(r'/\*(.*?)\*/', re.S)

    text = r.sub('', text)
    f = open('qhull.c', 'w')
    f.write(text)
    f.close()
finally:
    shutil.rmtree(tmp_dir)
