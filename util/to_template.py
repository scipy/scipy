#!/usr/bin/python

# takes source or interface file .xxx and produces templated .xxx.src
#   file for use in extending functionality to multiple precisions.
# the source file is assumed to use double precision D
#   for the first name of subroutines.
#   this is replaced  with template rules for all precisions


# <...>  is the template All blocks in a source file with names that
#         contain '<..>' will be replicated according to the
#         rules in '<..>'.

# The number of comma-separeted words in '<..>' will determine the number of
#   replicates.
 
# '<..>' may have two different forms, named and short. For example,

#named:
#   <p=d,s,z,c> where anywhere inside a block '<p>' will be replaced with
#  'd', 's', 'z', and 'c' for each replicate of the block.

#short:
#  <d,s,z,c>, a short form of the named, useful when no <p> appears inside 
#  a block.

#  Note that all <..> forms in a block must have the same number of
#    comma-separated entries. 

import string,os,sys
if sys.version[:3]>='2.3':
    import re
else:
    import pre as re

reg = re.compile(r"\ssubroutine\s(.+)\(.*\)")

try:
    file = sys.argv[1]
except IndexError:
    fid = sys.stdin
    outfile = sys.stdout
else:
    fid = open(file,'r')
    newname = file + '.src'
    outfile = open(newname,'w')

allstr = fid.read()
newstr = allstr.lower()
subs = reg.findall(newstr)
if subs is None:
    print "Nothing to do..."
    sys.exit(1)
for sub in subs:
    newstr = newstr.replace(sub,'<s,d,c,z>'+sub[1:])
    
newstr = newstr.replace('double precision','<real, double precision, complex, double complex>')

outfile.write(newstr)
