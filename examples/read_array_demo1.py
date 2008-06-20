#=========================================================================
# NAME: read_array_demo1
#
# DESCRIPTION: Examples to read 2 columns from a multicolumn ascii text
# file, skipping the first line of header. First example reads into
# 2 separate arrays. Second example reads into a single array. Data are
# then plotted.
#
# Here is the format of the file test.txt:
# --------
# Some header to skip
# 1 2 3
# 2 4 6
# 3 6 9
# 4 8 12
#
# USAGE:
# python read_array_demo1.py
#
# PARAMETERS:
#
# DEPENDENCIES:
# matplotlib (pylab)
# test.txt
#
#
# AUTHOR: Simon J. Hook
# DATE  : 09/23/2005
#
# MODIFICATION HISTORY:
#
# COMMENT:
#
#============================================================================

from scipy import *
from scipy.io import read_array
from pylab import *

def main():

    # First example, read first and second column from ascii file. Skip first
    # line of header.
    # Note use of (1,-1) in lines to skip first line and then read to end of file
    # Note use of (0,) in columns to pick first column, since its a tuple need trailing comma
    x=read_array("test.txt",lines=(1,-1), columns=(0,))
    y=read_array("test.txt",lines=(1,-1), columns=(1,))

    #Second example, read the file into a single arry
    z=read_array("test.txt",lines=(1,-1), columns=(0,2))

    # Plot the data
    plot(x,y,'r--',z[:,0],z[:,1])
    show()

# The one and only main function
if __name__ == "__main__":
    main()
