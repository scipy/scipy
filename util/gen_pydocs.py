#!/usr/bin/env python
"""
Script for generating scipy documentation using hacked version of pydoc.
Requires working installation of scipy.
Usage:

  python get_pydocs.py

This will generate a tree of html files into the 'pydoc' directory.
Check the end of this file to enable generating docs for various scipy
modules.

Bugs:
  - Documentation for some of the modules may not be generated.
    Fix: import these modules explicitely before write_docs() call.
  - When xxx/__init__.py does this:
     from foo import foo
    then docs for xxx.foo are not generated.
    Fix: call write_docs(sys.modules['xxx.foo'])

Author: Pearu Peterson <pearu@cens.ioc.ee>
Created: April 14, 2004
"""

import os
import sys
import my_pydoc as pydoc
from scipy_base.ppimport import disable
import scipy
# need to import scipy before disable to workaround circular import problems
disable()

output_dir = 'pydoc'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
os.chdir(output_dir)


def write_docs(module,prefix=None,_cache=[],_prefix=[]):
    top_call = ''
    if not _prefix or prefix is not None:
        if prefix is not None:
            if _prefix:
                _prefix[0] = prefix
            else: _prefix.append(prefix)
        else:
            _prefix.append(module.__name__.split('.')[0])
        if not os.path.exists(_prefix[0]):
            os.makedirs(_prefix[0])
        top_call = os.path.abspath(os.getcwd())
        os.chdir(_prefix[0])
    nofdots = _prefix[0].count('.')
    name = module.__name__

    if name.count('.')<nofdots:
        return
    if '.'.join(name.split('.')[:nofdots+1]) != _prefix[0] or name in _cache:
        return
    _cache.append(module.__name__)

    pydoc.writedoc(module)

    for n in dir(module):
        a = getattr(module,n)
        if type(a) is type(module):
            write_docs(a)
        elif str(type(a))=="<type 'fortran'>":
            a.__doc__

    if top_call:
        _cache.sort()
        d = {}
        for n in _cache:
            s = n.split('.')
            if not d.has_key(s[0]):
                d[s[0]] = []
            d1 = d[s[0]]
            while s:
                if len(s)==1: break
                fl = 1
                for dd in d1:
                    if dd.has_key(s[1]):
                        d1 = dd[s[1]]
                        s = s[1:]
                        fl = 0
                        break
                if fl:
                    d1.append({s[1]:[]})
                    s = s[1:]
        make_index(d)
        os.chdir(top_call)
        _cache[:] = []
        _prefix[:] = []

def make_index(d,depth=0,parent=[],fo = None):
    if depth==0:
        fo = open('index.html','w')

    if not d:
        return
    keys = d.keys()
    assert len(keys)==1,`keys`
    key = keys[0]
    fo.write('%s- "%s":%s.html' \
             % ('  '*depth+(depth and ' ' or ''),key,'.'.join(parent+[key])))
    if has_subtree(d[key]):
        fo.write('\n\n')
        for d1 in d[key]:
            make_index(d1,depth=depth+1,parent=parent+[key],fo = fo)
    else:
        c = []
        for d1 in d[key]:
            k = d1.keys()[0]
            c.append('"%s":%s.html' % (k,'.'.join(parent+[key,k])))
        if c:
            fo.write(' -> %s\n\n' % (', '.join(c)))
        else:
            fo.write('\n\n')

    if depth==0:
        fo.close()

def has_subtree(l):
    r = 0
    for d in l:
        if d[d.keys()[0]]:
            r = 1
            break
    return r

if 0:
    write_docs(scipy)

if 0:
    import scipy_base
    write_docs(scipy_base)

if 0:
    import scipy_test
    write_docs(scipy_test)

if 0:
    import scipy_distutils
    import scipy_distutils.system_info
    import scipy_distutils.fcompiler
    scipy_distutils.fcompiler.show_fcompilers()
    import scipy_distutils.mingw32ccompiler
    import scipy_distutils.lib2def
    import scipy_distutils.extension
    import scipy_distutils.core
    write_docs(scipy_distutils)

if 0:
    import weave
    write_docs(weave)

if 0:
    import gui_thread
    import gui_thread.wxPython_thread
    import gui_thread.wxBackgroundApp
    import gui_thread.main
    write_docs(gui_thread)
    write_docs(sys.modules['gui_thread.wxPython_thread'])
