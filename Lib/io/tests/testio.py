#!/usr/bin/env python

# This python script tests the numpyio module.
# also check out numpyio.fread.__doc__ and other method docstrings.

import scipy.numeric as Numeric
import numpyio
from MLab import rand
import os

# Generate some data
a = 255*rand(20)

# Open a file
fid = open("testfile.dat","w")

# Write the data as shorts
numpyio.fwrite(fid,20,a,Numeric.Int16)
fid.close()

# Reopen the file and read in data
fid = open("testfile.dat","r")
print "Don't worry about a warning regarding the number of bytes read.\n"
b = numpyio.fread(fid,1000000,Numeric.Int16,Numeric.Int)
fid.close()

if (Numeric.product(Numeric.equal(b,a.astype(Numeric.Int16))) == 0):
    print "Error in reading or writing..."

os.remove("testfile.dat")
