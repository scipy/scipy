#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

import tempfile
import subprocess
import os
import sys
import re
import shutil

from mako.template import Template

f = open('interpnd.pyx.in', 'r')
template = f.read()
f.close()

# Run templating engine
f = open('interpnd.pyx', 'w')
f.write(Template(template).render())
f.close()

tmp_dir = tempfile.mkdtemp()
try:
    # Run Cython
    dst_fn = os.path.join(tmp_dir, 'interpnd.c')
    ret = subprocess.call(['cython', '-I', '../..', '-o', dst_fn, 'interpnd.pyx'])
    if ret != 0:
        sys.exit(ret)

    # Strip comments
    f = open(dst_fn, 'r')
    text = f.read()
    f.close()

    r = re.compile(r'/\*(.*?)\*/', re.S)

    text = r.sub('', text)
    f = open('interpnd.c', 'w')
    f.write(text)
    f.close()
finally:
    shutil.rmtree(tmp_dir)
