#!/usr/bin/env python
import tempfile
import subprocess
import os
import sys
import re
import shutil

from mako.template import Template

f = open('interpnd0.pyx', 'r') # open unprocessed file
template = f.read()
f.close()
fn = 'interpnd.pyx' # processed file 

try:
    # Run templating engine
    f = open(fn, 'w')
    f.write(Template(template).render())
    f.close()

    # Run Cython
    dst_fn = 'interpnd.c'
    ret = subprocess.call(['cython', '-I', '../..', '-o', dst_fn, fn])
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
	os.unlink(fn)
