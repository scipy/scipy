#! /usr/bin/env python 
#

import sys
sys.path.append('..')

import cephes,Test

import glob, re

if sys.argv[-1] == '-b': reference=(0==0)
else: reference=(0==1)

for filenam in glob.glob('fncs_*.dat'):
    file=open(filenam,'r')
    _last=0
    while 1:
        testline=""
        while 1:
            line=file.readline()
            if line=="":
                _last=1
                break
            if line[-1:] == '\n': line = line[:-1]
            if line[-1:] == '\r': line = line[:-1]
            line=re.sub("#.*","",line)
            if re.search('\S',line):
                testline=testline+line
            else: break
        if testline != "":
            a=eval('Test.Test('+testline+',ref=reference)')
            rel_err=a.test()
            print a.name,
            if rel_err>0: print 'MAX_ERR: ', rel_err*100,'%'
            else: print 'PASS'
        if _last: break
    
