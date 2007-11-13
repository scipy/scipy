#!/usr/bin/env python

"""generate cephes_doc.h from included_functions.html"""


def parse(infile):
    d={}
    key=None
    val=''
    prev_line = ''
    for line in infile.readlines():
        if not line.strip():
            continue
        if line[0]=='<':
            if key and val:
                d[key]=val.strip()
                key,val=None,None
            if line[:4]=='<DT>':
                tok=line.split()
                tok=tok[-1].split('(')
                key=tok[0]
            elif line[:4]=='<DD>' and key:
                prev_line = prev_line[4:]
                tok = prev_line.split(' = ')
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
    ofile.write('#ifndef CEPHES_DOC_H\n')
    ofile.write('#define CEPHES_DOC_H\n')
    for key in keys:
        ofile.write('#define %s_doc "%s"\n'%(key,repr(d[key])[1:-1]))
    ofile.write('#endif /* CEPHES_DOC_H */\n')
    ofile.close()
