#!/usr/bin/env python

"""generate cephes_doc.h from included_functions.html"""

import string

def parse(infile):
    d={}
    key=None
    val=''
    prev_line = ''
    for line in infile.readlines():
        if not string.strip(line):
            continue
        if line[0]=='<':
            if key and val:
                d[key]=string.strip(val)
                key,val=None,None
            if line[:4]=='<DT>':
                tok=string.split(line)
                tok=string.split(tok[-1],'(')
                key=tok[0]
            elif line[:4]=='<DD>' and key:
                prev_line = prev_line[4:]
                tok = string.split(prev_line,' = ')
                val=tok[0]+'='+line[4:]
        else:
            if val:
                val=val+line
        prev_line = line

    return d

if __name__=="__main__":
    d = parse(open("docs/included_functions.html",'r'))
    keys = d.keys()
    keys.sort()
    ofile=open("cephes_doc.h",'w')
    for key in keys:
        ofile.write('#define %s_doc "%s"\n'%(key,repr(d[key])[1:-1]))
    ofile.close()
