from Numeric import *
from Scientific.IO.NetCDF import *
import time


def getUserName():
    try:
        import os, pwd, string
    except ImportError:
        return 'unknown user'
    pwd_entry = pwd.getpwuid(os.getuid())
    name = string.strip(string.splitfields(pwd_entry[4], ',')[0])
    if name == '':
        name = pwd_entry[0]
    return name

#
# Creating a file
#
file = NetCDFFile('test.nc', 'w', 'Created ' + time.ctime(time.time())
                  + ' by ' + getUserName())

file.title = "Just some useless junk"
file.version = 42

file.createDimension('xyz', 3)
file.createDimension('n', 20)
file.createDimension('t', None) # unlimited dimension

foo = file.createVariable('foo', Float, ('n', 'xyz'))
foo[:,:] = 0.
foo[0,:] = [42., 42., 42.]
foo[:,1] = 1.
foo.units = "arbitrary"
print foo[0]
print foo.dimensions

bar = file.createVariable('bar', Int, ('t', 'n'))
for i in range(10):
    bar[i] = i
print bar.shape

print file.dimensions
print file.variables

file.close()

#
# Reading a file
#
file = NetCDFFile('test.nc', 'r')

print file.dimensions
print file.variables

foo = file.variables['foo']
foo_array = foo[:]
foo_units = foo.units
print foo[0]

file.close()
